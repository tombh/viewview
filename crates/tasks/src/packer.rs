//! Find a not-terrible packing of Total Viewshed tiles over the entire globe.
//!
//! There are only 3 critical requirements:
//!   1  That each tile is perfectly square.
//!   2. That all tiles overlap. Or in other words, that there are no gaps, on land at least.
//!   3. That each tile is big enough to contain the theoretically longest line of sight of its
//!      highest point.

use std::{collections::VecDeque, panic};

use color_eyre::Result;
use geo::{BooleanOps as _, Contains as _, GeodesicArea as _, Intersects as _, Within as _};
use rstar::PointDistance as _;

use crate::projector::LatLonCoord;

/// The minimum elevation inside a tile to start considering larger tiles from. Another way of
/// looking at this is rather the minimum _width_ of tile, as tile widths are primarilly dictated
/// by the maximum point of elevation inside them. We don't use 0m because in fact the lowest
/// elevations can be below zero, which, although geographically valid, cannot be used to calculate
/// a highest theoretical maximum line of sight. So if negative and zero values aren't useful, then
/// we just choose a reasonable low elevation that generates tile sizes worthy of going to the
/// effort of calculating for.
const MINIMUM_HEIGHT: i32 = 100;

/// The maximum radius in meters of the moving window view.
///
/// Let's assume that the longest line of sight is around 600km. Then our window radius needs to be
/// at least that. But then if we're analysing a point on the edge of the window it needs to also
/// see points 600km away. So the minimum window size is then 2*600km. So let's add yet another
/// 600km just to be safe. It still runs nice and fast.
const WINDOW_RADIUS: f64 = 1_800_000.0;
/// The distance that each window moves.
const WINDOW_STEP: f64 = WINDOW_RADIUS / 2.0f64;

/// A single point of elevation that represents the highest elevation within the resolution range of
/// the point. The resolution is defined by another process, `max_subtile.rs`.
type ElevationRstar = rstar::primitives::GeomWithData<geo::Coord, i32>;
/// How we store tiles for fast lookups.
type TileRstar = rstar::primitives::GeomWithData<geo::Coord, crate::tile::Tile>;

/// A packer of tiles.
pub struct Packer {
    /// Config from the CLI.
    config: crate::config::Packer,
    /// For projecting between different coordinate systems. There are 3:
    ///   1. Conventional lon/lat.
    ///   2. AEQD anchored to the centre of the current packing "window". An expensive part of the
    ///      packer is finding all the elevations within a tile. Ideally we would reproject all
    ///      the millions of elevations to be anchored on the tile, but I think as long as the
    ///      current window is within a reasonable distance, say 2000km, then there isn't a
    ///      significant loss of accuracy.
    ///      
    ///      TODO:
    ///      First, find out what the error margin is. And then if it is a problem, explore
    ///      using a lon/lat AABB to find a reduced set of points that can be reprojected to the
    ///      tile's centre.
    ///   3. AEQD anchored to the cnetre of a given tile. This is necessary for things like
    ///      constructing the tile's actual corners. The simple act of adding a distance to a point
    ///      must be done in as local as possible projection.
    projector: crate::projector::Convert,
    /// An unchanging canonical reference of all the world's maximum elevation points.
    canonical: rstar::RTree<ElevationRstar>,
    /// The window is a view onto a more manageable subset of the canonical data. It is mainly for
    /// faster lookups but also helps with projection accuracy.
    window: rstar::RTree<ElevationRstar>,
    /// A stack of points, each of which must at some point be proven to fall within a tile.
    stack: VecDeque<ElevationRstar>,
    /// All the currently found tiles for the _window_.
    tiles: rstar::RTree<TileRstar>,
    /// All the currently found tiles for the _world_.
    all: rstar::RTree<TileRstar>,
}

impl Packer {
    /// Instantiate.
    pub fn new(config: crate::config::Packer) -> Result<Self> {
        Ok(Self {
            config,
            canonical: Self::build_canonical_rtree()?,
            window: rstar::RTree::new(),
            stack: VecDeque::new(),
            tiles: rstar::RTree::new(),
            all: rstar::RTree::new(),
            projector: crate::projector::Convert {
                base: crate::projector::LatLonCoord(geo::Coord::zero()),
            },
        })
    }

    /// Run for the entire globe.
    pub fn run_all(&mut self) -> Result<()> {
        let default_start = (-180.0f64, 90.0f64);
        let mut centre = LatLonCoord(self.config.start.unwrap_or(default_start).into());

        let mut steps = 0u32;
        while (-90.0f64..=90.0f64).contains(&centre.0.y) {
            if let Some(stop_at_step) = self.config.steps
                && steps >= stop_at_step
            {
                break;
            }

            let longtitude_step = Self::calculate_longtitude_step(centre.0.y);
            let mut is_last_longitude = false;
            while (-180.0f64..=180.0f64).contains(&centre.0.x) {
                if let Some(stop_at_step) = self.config.steps
                    && steps >= stop_at_step
                {
                    break;
                }

                tracing::debug!("Calculating tiles for window centred at: {centre:?}");

                self.run_one(centre)?;
                self.nudge_extend_window_tiles()?;
                self.move_tiles_to_all()?;
                self.remove_inefficient_tiles()?;
                self.save_tiles_as_geojson()?;
                if is_last_longitude {
                    break;
                }

                centre = LatLonCoord(geo::coord! {
                    x: centre.0.x + longtitude_step,
                    y: centre.0.y
                });

                if centre.0.x > 180.0f64 {
                    centre.0.x = 180.0f64;
                    is_last_longitude = true;
                }

                steps += 1;
            }
            centre = Self::start_new_latitude(centre)?;
        }

        Ok(())
    }

    /// Move tiles generated in a window to the global store.
    fn move_tiles_to_all(&mut self) -> Result<()> {
        for tile in &self.tiles {
            let tile_lonlat = TileRstar::new(self.projector.to_degrees(*tile.geom())?.0, tile.data);
            self.all.insert(tile_lonlat);
        }
        self.tiles = rstar::RTree::new();
        Ok(())
    }

    /// Calculate the number of degrees to go Eastward to reach the next window.
    fn calculate_longtitude_step(latitude: f64) -> f64 {
        let degrees_per_meter = 1.0 / crate::projector::Convert::meters_per_degree(latitude);
        f64::from(degrees_per_meter) * WINDOW_STEP
    }

    /// Do a sort of "carriage return" to the next latitude South.
    fn start_new_latitude(current: LatLonCoord) -> Result<LatLonCoord> {
        let projector = crate::projector::Convert { base: current };
        let down = projector.to_degrees(geo::coord! {
            x: 0.0f64,
            y: -WINDOW_STEP,
        })?;
        Ok(LatLonCoord((-180.0, down.0.y).into()))
    }

    /// Run a single step.
    pub fn run_one(&mut self, centre: LatLonCoord) -> Result<()> {
        self.build_window(centre);
        self.iterate()?;

        if self.config.one.is_some() {
            self.nudge_extend_window_tiles()?;
            self.remove_inefficient_tiles()?;
            self.save_tiles_as_geojson()?;
        }

        Ok(())
    }

    /// Build a view onto a subset of the total data. This primarily allows for faster lookups of
    /// relevant data, and also hopefully adds a bit to the accuracy as all points in the window
    /// are projected to an AEQD projection anchored on the centre of the window.
    fn build_window(&mut self, centre: LatLonCoord) {
        let mut window = Vec::new();
        self.stack = VecDeque::new();
        self.projector = crate::projector::Convert { base: centre };

        tracing::debug!("Ordering around coordinate: {centre:?}");
        for canonical_point in self.canonical.nearest_neighbor_iter(&centre.0) {
            let lonlat = crate::projector::LatLonCoord(*canonical_point.geom());
            #[expect(clippy::panic, reason = "This should be considered a bug.")]
            let projected = self
                .projector
                .to_meters(lonlat)
                .unwrap_or_else(|_| panic!("Couldn't project point: {lonlat:?}"));
            let height = canonical_point.data;
            let point = ElevationRstar::new(projected, height);
            window.push(point);

            let distance = point.geom().distance_2(&geo::Coord::zero()).sqrt();
            if distance < WINDOW_STEP {
                self.stack.push_back(point);
            }
            if distance > WINDOW_RADIUS {
                break;
            }
        }

        // TODO: remember to also reproject the tile RTree when moving the window.

        self.window = rstar::RTree::bulk_load(window);
    }

    /// The centre of the current window.
    const fn current_window_centre(&self) -> LatLonCoord {
        self.projector.base
    }

    /// Load the acceleration structire of maximum elevation subtiles for the whole world.
    fn build_canonical_rtree() -> Result<rstar::RTree<ElevationRstar>> {
        let mut points: Vec<ElevationRstar> = Vec::new();
        // TODO: make the file name configurable.
        let subtiles = crate::max_subtile::Subtiler::load("max_subtiles.bin")?;
        let total_points = subtiles.len();

        tracing::info!("Loading {} max subtiles.", total_points);
        for subtile in subtiles {
            // TODO: should this be using the LonLatCoord type?
            points.push(ElevationRstar::new(
                geo::coord! { x: subtile.lon.into(), y: subtile.lat.into()},
                subtile.max_height,
            ));
        }

        tracing::info!(
            "Creating canonical RTree from {} max elevation points",
            total_points
        );
        let tree: rstar::RTree<_> = rstar::RTree::bulk_load(points);

        Ok(tree)
    }

    /// Loop over all the points in a window.
    fn iterate(&mut self) -> Result<()> {
        let original_points: Vec<ElevationRstar> = self.stack.clone().into();
        let mut zeroes = Vec::new();
        while let Some(point) = self.get_next_point() {
            if point.data <= MINIMUM_HEIGHT {
                zeroes.push(point);
                continue;
            }
            self.find_non_overlapping_tile_for_point(point)?;
        }

        self.stack = zeroes.into();
        while let Some(point) = self.get_next_point() {
            self.find_non_overlapping_tile_for_point(point)?;
        }

        self.iterate_underlappers(original_points)?;

        Ok(())
    }

    /// Search for the best non-overlapping tile for the given point.
    fn find_non_overlapping_tile_for_point(&mut self, centre: ElevationRstar) -> Result<()> {
        tracing::debug!("Finding non-overlapping tile for point: {centre:?}");
        let mut highest = centre;
        loop {
            let tile = self.find_minimum_tile_for_point(&centre, highest.data)?;
            let Some((candidate, covered_points)) = self.find_highest_point_in_tile(&tile)? else {
                // Ignore tiles with just zeroes in them.
                return Ok(());
            };

            if candidate.data > highest.data {
                tracing::trace!(
                    "Found new highest elevation: {} > {}",
                    candidate.data,
                    highest.data
                );
                highest = candidate;
            } else {
                let overlap_amount = self.calculate_overlap_amount(&tile)?;
                if overlap_amount == 0.0 {
                    self.add_tile(tile, &covered_points)?;
                    return Ok(());
                }

                tracing::trace!("Overlap: {overlap_amount}");
                return Ok(());
            }
        }
    }

    /// Search for the best tile for the given point regardless of other tiles.
    fn find_any_tile_for_point(
        &self,
        centre: ElevationRstar,
        start_with_highest: Option<i32>,
    ) -> Result<Option<(crate::tile::Tile, Vec<ElevationRstar>)>> {
        tracing::debug!("Finding any tile for point: {centre:?}");
        let mut highest = start_with_highest.unwrap_or(centre.data);
        loop {
            let tile = self.find_minimum_tile_for_point(&centre, highest)?;
            let Some((candidate, covered_points)) = self.find_highest_point_in_tile(&tile)? else {
                // Ignore tiles with just zeroes in them.
                return Ok(None);
            };

            if candidate.data > highest {
                tracing::trace!(
                    "Found new highest elevation: {} > {}",
                    candidate.data,
                    highest
                );
                highest = candidate.data;
            } else {
                return Ok(Some((tile, covered_points)));
            }
        }
    }

    /// Find all the points in a tile.
    //
    // TODO: Investigate doing this all in lon/lat. It may save some reprojection.
    fn find_points_in_tile(
        &self,
        tile: &crate::tile::Tile,
    ) -> Result<Vec<rstar::primitives::GeomWithData<geo::Coord, i32>>> {
        let tile_polygon = tile.to_polygon_metric(self.current_window_centre())?;
        let mut covered_points: Vec<ElevationRstar> = self
            .window
            .locate_in_envelope_intersecting(&tile.to_aabb_metric(self.current_window_centre())?)
            .copied()
            .collect();
        covered_points.retain(|point| point.geom().is_within(&tile_polygon));
        Ok(covered_points)
    }

    /// Calculate how much a tile overlaps with its neighbours.
    fn calculate_overlap_amount(&self, tile: &crate::tile::Tile) -> Result<f64> {
        let candidate = tile.to_polygon_lonlat()?;
        let area = candidate.geodesic_area_unsigned();
        let mut overlap = 0.0f64;
        let mut searched = 0usize;
        for neighbour_basic in self.tiles.nearest_neighbor_iter(&tile.centre.0) {
            let neighbour = neighbour_basic.data.to_polygon_lonlat()?;
            let intersection = candidate.intersection(&neighbour);
            overlap += intersection.geodesic_area_unsigned();
            searched += 1;
            if searched > 1000 {
                break;
            }
        }
        Ok(overlap / area)
    }

    /// Find points that have not yet been associated with a tile.
    fn find_underlappers(
        &self,
        original_points: Vec<ElevationRstar>,
    ) -> Result<Vec<ElevationRstar>> {
        let mut underlappers = original_points;
        for tile in &self.tiles {
            let covered_points = self.find_points_in_tile(&tile.data)?;
            underlappers.retain(|point| !covered_points.contains(point));
        }

        Ok(underlappers)
    }

    /// Loop over all points that have not yet been covered by a tile.
    fn iterate_underlappers(&mut self, original_points: Vec<ElevationRstar>) -> Result<()> {
        let underlappers = self.find_underlappers(original_points)?;
        self.stack = underlappers.into();

        while let Some(point) = self.get_next_point() {
            self.handle_underlapper(point)?;
        }

        Ok(())
    }

    /// Handle a point that isn't yet contained within a tile.
    ///
    /// Make 2 tiles:
    ///   1. A tile for the point, even though it will overlap another tile.
    ///   2. Find the nearest tile to the underlapping point and extend it to cover the point.
    ///
    /// And add the one that introduces the least surface area.
    fn handle_underlapper(&mut self, point: ElevationRstar) -> Result<()> {
        let Some((own_tile, own_covered_points)) = self.find_any_tile_for_point(point, None)?
        else {
            tracing::trace!("No tile found");
            return Ok(());
        };
        let own_surface_area = own_tile.surface_area()?;

        if !self.attempt_increasing_nearest_tile(point, own_surface_area)? {
            self.add_tile(own_tile, &own_covered_points)?;
        }

        Ok(())
    }

    /// Increase the size of the nearest tile to the point such that it covers the point. If that
    /// size increase is less than making a dedicated tile for the point, than accept the increased
    /// size of the tile.
    fn attempt_increasing_nearest_tile(
        &mut self,
        point: ElevationRstar,
        surface_area_to_beat: f32,
    ) -> Result<bool> {
        if self.tiles.size() == 0 {
            return Ok(false);
        }
        let (mut nearest_tile, old_width) = self.engulf_tile_to_point(&point)?;
        let nearest_tile_starting_surface_area = nearest_tile.data.surface_area()?;
        tracing::trace!("🏝️ Nearest tile for underlapping point ({point:?}): {nearest_tile:?}");

        self.shrink_wrap_tile_to_point(&mut nearest_tile, point)?;
        let nearest_tile_covered_points = self.ensure_tile_is_big_enough(&mut nearest_tile)?;

        let nearest_tile_surface_increase =
            nearest_tile.data.surface_area()? - nearest_tile_starting_surface_area;
        if nearest_tile_surface_increase > surface_area_to_beat {
            return Ok(false);
        }

        tracing::trace!(
            "Changed nearest tile's width: {old_width} to {}",
            nearest_tile.data.width
        );
        assert!(
            nearest_tile.data.width > old_width,
            "Resized tile didn't increase in size, probably a bug."
        );

        self.tiles.pop_nearest_neighbor(point.geom());
        self.add_tile(nearest_tile.data, &nearest_tile_covered_points)?;
        Ok(true)
    }

    /// Increase a tile's width such that it contains the given point.
    ///
    /// Overshoots to ensure no false positives between projection changes.
    fn engulf_tile_to_point(&self, point: &ElevationRstar) -> Result<(TileRstar, f32)> {
        let overshoot = 1.01;
        let Some(nearest_tile_reference) = self.tiles.nearest_neighbor(point.geom()) else {
            color_eyre::eyre::bail!("No nearest tile found for point: {point:?}");
        };
        let mut nearest_tile = *nearest_tile_reference;
        let old_width = nearest_tile.data.width;

        let tile_projecter = crate::projector::Convert {
            base: nearest_tile.data.centre,
        };

        // Relative distances must be done in a projection anchored to the point where the relative
        // distance begins.
        let point_lonlat = self.projector.to_degrees(*point.geom())?;
        let point_anchored_to_tile = tile_projecter.to_meters(point_lonlat)?;

        let distance = point_anchored_to_tile.x.hypot(point_anchored_to_tile.y);
        tracing::trace!("Point is {distance}m from its nearest tile.");
        #[expect(
            clippy::cast_possible_truncation,
            clippy::as_conversions,
            reason = "I think this is the only way?"
        )]
        let first_width = (distance as f32 * 2.0) * overshoot;
        nearest_tile.data.width = first_width;

        Ok((nearest_tile, old_width))
    }

    /// Slowly shrink a tile until we find the last size within which `is_within()` returns true
    /// for the given point.
    fn shrink_wrap_tile_to_point(&self, tile: &mut TileRstar, point: ElevationRstar) -> Result<()> {
        let first_width = tile.data.width;
        let step = 100.0;
        loop {
            let polygon = tile.data.to_polygon_metric(self.current_window_centre())?;
            if !point.geom().is_within(&polygon) {
                #[expect(
                    clippy::float_cmp,
                    reason = "`first_width` is the original value, not the result of a calculation."
                )]
                {
                    assert!(
                        tile.data.width != first_width,
                        "Resized nearest tile never contained the underlapping point"
                    );
                };
                tile.data.width += step;
                break;
            }
            tile.data.width -= step;
        }

        tracing::trace!(
            "New width required for nearest tile to contain point: {}",
            tile.data.width
        );

        Ok(())
    }

    /// When a tile is resized after its initial creation, we need to check that the new size
    /// hasn't covered a new higher elevation that would cause the tile to need to be resized yet
    /// again.
    #[expect(
        clippy::panic,
        reason = "
          Panicking is better as this is just a one-off script, so its better to just fix the bug.
        "
    )]
    fn ensure_tile_is_big_enough(&self, tile: &mut TileRstar) -> Result<Vec<ElevationRstar>> {
        let Some((tile_highest, mut tile_covered_points)) =
            self.find_highest_point_in_tile(&tile.data)?
        else {
            panic!(
                "Should be impossible because tile was already created with non-zero elevations"
            );
        };

        let start_at_height = tile_highest.data;
        let Some((resized_tile, resized_covered_points)) =
            self.find_any_tile_for_point(ElevationRstar::new(*tile.geom(), start_at_height), None)?
        else {
            panic!("We can't recover from this. It's probably a bug");
        };
        if resized_tile.width > tile.data.width {
            tracing::debug!(
                "Resized tile was resized _again_ due to containing a higher highest point:"
            );
            tracing::debug!("  From: {}m to {}m", tile.data.width, resized_tile.width,);
            tile.data.width = resized_tile.width;
            tile_covered_points = resized_covered_points;
        } else {
            tracing::trace!(
                "Resized tile doesn't contain higher points, so doesn't need to be resized again."
            );
        }

        Ok(tile_covered_points)
    }

    /// Add a single tile to the list of all tiles that covered the world's land.
    fn add_tile(
        &mut self,
        tile: crate::tile::Tile,
        covered_points: &[ElevationRstar],
    ) -> Result<()> {
        self.tiles.insert(TileRstar::new(
            tile.centre_metric(self.current_window_centre())?,
            tile,
        ));
        tracing::info!("Added tile: {tile:?}");

        if self.config.one.is_some() {
            self.save_tiles_as_geojson()?;
        }

        let before = self.stack.len();
        self.stack.retain(|point| !covered_points.contains(point));
        let after = self.stack.len();
        tracing::trace!(
            "{} points in tile, {} removed, {} points remaining.",
            covered_points.len(),
            before - after,
            self.stack.len()
        );
        Ok(())
    }

    /// Pop the next point off the stack.
    fn get_next_point(&mut self) -> Option<ElevationRstar> {
        self.stack.pop_front()
    }

    /// Find the minimum bounding box that would allow calculating all the lines of site based on
    /// how far the highest point can see.
    fn find_minimum_tile_for_point(
        &self,
        point: &ElevationRstar,
        highest: i32,
    ) -> Result<crate::tile::Tile> {
        let width = Self::minimum_tile_size_for_elevation(highest);
        tracing::trace!("Minimum extent for {point:?}: {width}");
        let centre = self.projector.to_degrees(*point.geom())?;
        let tile = crate::tile::Tile { centre, width };
        Ok(tile)
    }

    /// Given a height, how far is the longest theoretical line of sight from it.
    // TODO:
    //   Add a subtile of meters calculated at the furthest point?
    //   I don't think so, I can't think of an intuitive explanation though.
    fn minimum_tile_size_for_elevation(elevation_i32: i32) -> f32 {
        if elevation_i32 <= MINIMUM_HEIGHT {
            // Strike a balance between lots of tiny tiles or fewer bigger tiles.
            return 36_000.0;
        }

        #[expect(
            clippy::as_conversions,
            clippy::cast_precision_loss,
            reason = "Is there another way?"
        )]
        let elevation = elevation_i32 as f32;

        (2.0 * crate::projector::EARTH_RADIUS * 1000.0)
            .mul_add(elevation, elevation.powi(2))
            .sqrt()
            * 2.0
    }

    /// Find the highest DEM point in a candidate tile.
    fn find_highest_point_in_tile(
        &self,
        tile: &crate::tile::Tile,
    ) -> Result<Option<(ElevationRstar, Vec<ElevationRstar>)>> {
        let start = std::time::Instant::now();
        let mut highest = None;
        let mut covered_points = Vec::new();

        // The AABB is quicker to search within but is also slightly rotated in AEQD projections
        // so we need a second step to clip points outide of the tile's real extent.
        let tile_aabb = tile.to_aabb_metric(self.current_window_centre())?;
        let aabb_points = self.window.locate_in_envelope_intersecting(&tile_aabb);

        let covered_points_polygon = tile.to_polygon_metric(self.current_window_centre())?;
        let highest_points_polygon = tile.to_polygon_metric(self.current_window_centre())?;

        let mut current_highest = i32::MIN;
        let mut is_non_zero_detected = false;
        for point in aabb_points {
            if !point.geom().is_within(&highest_points_polygon) {
                continue;
            }
            if point.data != 0i32 {
                is_non_zero_detected = true;
            }
            if point.data > current_highest {
                current_highest = point.data;
                highest = Some(*point);
            }
            if point.geom().is_within(&covered_points_polygon) {
                covered_points.push(*point);
            }
        }

        if !is_non_zero_detected {
            return Ok(None);
        }

        if current_highest > i32::MIN {
            tracing::trace!(
                "Found highest point {:?} in {}μs",
                current_highest,
                start.elapsed().as_micros(),
            );
        }

        let Some(some_highest) = highest else {
            #[expect(clippy::panic, reason = "There's no way to recover")]
            {
                panic!("Tile candidate doesn't contain any data?");
            }
        };

        Ok(Some((some_highest, covered_points)))
    }

    /// Remove any tiles that are completely covered by other tiles.
    fn remove_inefficient_tiles(&mut self) -> Result<()> {
        let mut discarders = Vec::new();
        let mut keepers = Vec::new();
        let source = match self.config.one {
            Some(_) => &self.tiles,
            None => &self.all,
        };
        tracing::info!("Removing nested tiles from {} tiles...", source.size());

        for tile in source {
            if discarders.contains(tile) {
                continue;
            }

            let distance = tile.data.distance_from(self.current_window_centre())?;
            if distance > WINDOW_RADIUS * 1.5f64 {
                keepers.push(*tile);
                continue;
            }

            if self.is_tile_nested(tile, &discarders)? {
                discarders.push(*tile);
            } else {
                keepers.push(*tile);
            }
        }

        let count = match self.config.one {
            Some(_) => {
                self.tiles = rstar::RTree::bulk_load(keepers.clone());
                self.tiles.size() - keepers.len()
            }
            None => {
                self.all = rstar::RTree::bulk_load(keepers.clone());
                self.all.size() - keepers.len()
            }
        };

        tracing::info!("Removed {count} nested tiles.");
        Ok(())
    }

    /// Is the tile completely inside another tile?
    fn is_tile_nested(&self, tile: &TileRstar, discarders: &[TileRstar]) -> Result<bool> {
        let source = match self.config.one {
            Some(_) => &self.tiles,
            None => &self.all,
        };
        let polygon = tile.data.to_polygon_lonlat()?;
        let mut super_tile = geo::MultiPolygon::empty();
        let iterator = source.nearest_neighbor_iter(tile.geom()).enumerate();
        for (count, neighbour) in iterator {
            if neighbour == tile {
                continue;
            }
            if discarders.contains(neighbour) {
                continue;
            }
            if count > 50 {
                break;
            }

            let neighbour_polygon = neighbour.data.to_polygon_lonlat()?;
            if neighbour_polygon.intersects(&polygon) {
                super_tile = super_tile.union(&neighbour_polygon);
            }
        }

        Ok(super_tile.contains(&polygon))
    }

    /// Where the source of completed tiles are.
    const fn source_of_tiles(
        &self,
    ) -> &rstar::RTree<rstar::primitives::GeomWithData<geo::Coord, crate::tile::Tile>> {
        match self.config.one {
            Some(_) => &self.tiles,
            None => &self.all,
        }
    }

    /// Extend each tile by the resolution of the underlying max subtile data.
    ///
    /// Consider the case where 2 tiles perfectly align in the sense that one starts at the next
    /// row or column of points where the adjacent tile finished. This would mean that there is gap
    /// of exactly one max-subtile-resolution unit between the 2 tiles edges. This gap is _degree_
    /// based so the actual metric distance varies by latitude. And couple that with the fact that
    /// tiles have to be perfectly square, we can't just increase problematic tile edges by a
    /// constant amount.
    fn nudge_extend_window_tiles(&mut self) -> Result<()> {
        tracing::info!(
            "Nudge extending tiles ({}) in window to ensure they overlap.",
            self.tiles.size()
        );

        // This depends on the resolution set in the max-subtile pre-process step. It's not good
        // that we hardcode it here :/
        let subtile_resolution = 10.0;
        // Just add a little more for safe measure. Thought the better aproach would be to take the
        // latitude of the tile corner that is furthest from the equator.
        let magic = 1.5;

        let mut tiles = Vec::new();
        for mut tile in self.tiles.iter().copied() {
            let latitude = tile.data.centre.0.y;
            let extension = (crate::projector::Convert::meters_per_degree(latitude)
                / subtile_resolution)
                * magic;
            tile.data.width += extension;
            tiles.push(tile);
            self.ensure_tile_is_big_enough(&mut tile)?;
        }

        self.tiles = rstar::RTree::bulk_load(tiles);

        Ok(())
    }

    /// Save the tiles to [`GeoJSON`].
    fn save_tiles_as_geojson(&self) -> Result<()> {
        let mut features: Vec<geo::Geometry> = Vec::new();
        let path = "static/tiles.json";
        for tile in self.source_of_tiles() {
            let geojson_tile = tile.data.to_polygon_lonlat()?;
            features.push(geojson_tile.into());
        }

        tracing::info!("Saving {} tiles to: {path}", self.source_of_tiles().size());
        let feature_collection = geojson::FeatureCollection::from(&features.into());
        let json = geojson::GeoJson::from(feature_collection);
        std::fs::write(path, json.to_string())?;

        Ok(())
    }

    /// Save tiles and points [`GeoJSON`].
    #[expect(dead_code, reason = "Useful for debugging.")]
    fn save_tiles_and_points(
        maybe_tiles: Option<Vec<crate::tile::Tile>>,
        points: &[LatLonCoord],
    ) -> Result<()> {
        let mut features: Vec<geo::Geometry> = Vec::new();

        let path = "static/tiles.json";

        if let Some(tiles) = maybe_tiles.clone() {
            for tile in tiles {
                let geojson_tile = tile.to_polygon_lonlat()?;
                features.push(geojson_tile.into());
            }
        }

        let mut geojson_points = Vec::new();
        for point in points {
            geojson_points.push(geo::Point::from(point.0));
        }
        let multipoint = geo::MultiPoint::new(geojson_points);
        features.push(multipoint.into());

        let feature_collection = geojson::FeatureCollection::from(&features.into());

        tracing::info!(
            "DEBUG: saving {} tiles and {} points to: {path}",
            maybe_tiles.unwrap_or_default().len(),
            points.len()
        );
        let json = geojson::GeoJson::from(feature_collection);
        std::fs::write(path, json.to_string())?;

        Ok(())
    }

    /// Save polygons for debugging.
    #[expect(dead_code, reason = "Useful for debugging.")]
    fn save_polygons(polygons: &[geo::MultiPolygon]) -> Result<()> {
        let path = "static/tiles.json";
        let mut features: Vec<geo::Geometry> = Vec::new();
        for polygon in polygons.iter().cloned() {
            features.push(polygon.into());
        }

        let feature_collection = geojson::FeatureCollection::from(&features.into());

        tracing::info!("DEBUG: saving {} polygons to: {path}", polygons.len());
        let json = geojson::GeoJson::from(feature_collection);
        std::fs::write(path, json.to_string())?;

        Ok(())
    }
}
