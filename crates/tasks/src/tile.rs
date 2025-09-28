//! A single computable tile for the Total Viewshed algorithm.

use color_eyre::{Result, eyre::ContextCompat as _};
use geo::{BoundingRect as _, GeodesicArea as _, polygon};
use rstar::PointDistance as _;

use crate::projector::LatLonCoord;

/// The tile data itself.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tile {
    /// The centre of the tile.
    pub centre: crate::projector::LatLonCoord,
    /// The width of the tile. Therefore, not the distance from the centre to an edge.
    pub width: f32,
}

impl Tile {
    /// Convert the width to the shortest distance between the centre and an edge.
    fn offset(&self) -> f64 {
        f64::from(self.width / 2.0)
    }

    /// The centre coordinate reprojected to the given metric projection.
    pub fn centre_metric(&self, anchor: crate::projector::LatLonCoord) -> Result<geo::Coord> {
        let projecter = crate::projector::Convert { base: anchor };
        projecter.to_meters(self.centre)
    }

    /// The bounding box for the tile in lon/lat coordinates. Note that this is not a bounding box
    /// in the sense that it preserves a true metric square. It is "square" within degree-based
    /// units.
    fn to_bbox_lonlat(
        self,
        fudge: f32,
    ) -> Result<(
        crate::projector::LatLonCoord,
        crate::projector::LatLonCoord,
        crate::projector::LatLonCoord,
        crate::projector::LatLonCoord,
    )> {
        let projecter = crate::projector::Convert { base: self.centre };
        let offset = self.offset() * f64::from(1.0 + fudge);
        let top_left = projecter.to_degrees(geo::coord! { x: -offset, y: offset })?;
        let top_right = projecter.to_degrees(geo::coord! { x: offset, y: offset })?;
        let bottom_right = projecter.to_degrees(geo::coord! { x: offset, y: -offset })?;
        let bottom_left = projecter.to_degrees(geo::coord! { x: -offset, y: -offset })?;

        Ok((top_left, top_right, bottom_right, bottom_left))
    }

    /// The bounding box for the tile in metric coordinates. Note that the bouding box is actually
    /// just exactly the same as the tile itself. Perhaps there's a better name to use?
    fn to_bbox_metric(
        self,
        anchor: crate::projector::LatLonCoord,
        fudge: f32,
    ) -> Result<(geo::Coord, geo::Coord, geo::Coord, geo::Coord)> {
        let bbox_lonlat = self.to_bbox_lonlat(fudge)?;

        let anchor_projector = crate::projector::Convert { base: anchor };
        Ok((
            anchor_projector.to_meters(bbox_lonlat.0)?,
            anchor_projector.to_meters(bbox_lonlat.1)?,
            anchor_projector.to_meters(bbox_lonlat.2)?,
            anchor_projector.to_meters(bbox_lonlat.3)?,
        ))
    }

    /// The Axis-Aligned Bounding Box for the tile. Note that this has an unintuitive shape. The
    /// metric projection's 0,0 coord is still at a "polar" point, therefore axis-alignment rarely
    /// follows the axes of the tile. Nevertheless the AABB is still useful for quicker first pass
    /// lookups of containing points. A follow up `is_within()` for each found point can then be
    /// used to get the exact contents.
    pub fn to_aabb_metric(
        self,
        anchor: crate::projector::LatLonCoord,
        fudge: f32,
    ) -> Result<rstar::AABB<geo::Coord>> {
        let bbox = self
            .to_polygon_metric(anchor, fudge)?
            .bounding_rect()
            .context(format!("Couldn't find bbox for tile: {self:?}"))?;
        let aabb = rstar::AABB::from_corners(bbox.min(), bbox.max());

        Ok(aabb)
    }

    /// Lon/lat polygon of tile.
    pub fn to_polygon_lonlat(self) -> Result<geo::Polygon> {
        let bbox = self.to_bbox_lonlat(0.0)?;
        Ok(polygon![bbox.0.0, bbox.1.0, bbox.2.0, bbox.3.0, bbox.0.0])
    }

    /// Hack to get around tiles that don't overlap.
    pub fn to_polygon_lonlat_fudged(self, fudge: f32) -> Result<geo::Polygon> {
        let bbox = self.to_bbox_lonlat(fudge)?;
        Ok(polygon![bbox.0.0, bbox.1.0, bbox.2.0, bbox.3.0, bbox.0.0,])
    }

    /// Metric polygon of tile.
    pub fn to_polygon_metric(
        self,
        anchor: crate::projector::LatLonCoord,
        fudge: f32,
    ) -> Result<geo::Polygon> {
        let bbox = self.to_bbox_metric(anchor, fudge)?;
        Ok(polygon![bbox.0, bbox.1, bbox.2, bbox.3, bbox.0])
    }

    /// The surface area covered by the tile.
    pub fn surface_area(self) -> Result<f32> {
        let polygon = self.to_polygon_lonlat()?;
        #[expect(
            clippy::as_conversions,
            clippy::cast_possible_truncation,
            reason = "Is there another way?"
        )]
        Ok(polygon.geodesic_area_unsigned() as f32)
    }

    /// Calculate the distance in meters of the tile from the given point.
    pub fn distance_from(&self, point_lonlat: LatLonCoord) -> Result<f64> {
        let projector = crate::projector::Convert { base: point_lonlat };
        let point = projector.to_meters(self.centre)?;

        Ok(self.centre_metric(point_lonlat)?.distance_2(&point).sqrt())
    }
}
