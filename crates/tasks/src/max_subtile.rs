//! These subtiles are a kind of acceleration structure for being able to more quickly find the
//! highest points on the planet.
//!
//! We take a normal DEM tile, split it into a certain grid, say 10x10, and find the maximum
//! elevation in each tile of that 10x10 grid. We then record the lon/lat centre of each grid tile.
//!
//! This should be run on a large number of DEM tiles, typically in fact the whole world.

use std::io::Write as _;

use color_eyre::Result;

/// The number of arc seconds in a degree.
const ARCSEC_PER_DEG: f32 = 3600.0;

/// A `MaxSubTile` is sub tile of a DEM tile, that contains nothing other than the maximum height
/// of the subtile.
#[repr(C)]
#[derive(Copy, Clone, Debug, PartialEq, bytemuck::Zeroable, bytemuck::Pod)]
pub struct MaxSubTile {
    /// Longtitude
    lon: f32,
    /// Latitude
    lat: f32,
    /// The maximum height within the subtile region.
    max_height: i32,
}

pub struct Subtiler {
    /// Keep track of _all_ the subtiles on the planet.
    subtiles: Vec<MaxSubTile>,
    /// The width of the underlying DEM tile.
    tile_width: usize,
    /// The width of each subtile. Though note that the final subtiles of each row and colum may be
    /// less than this.
    subtile_width: usize,
    /// Keep track of how many invalid subtiles were found.
    pub invalid_count: u64,
}

impl Subtiler {
    /// Instantitate.
    pub fn new(factor: u32, width: usize) -> Result<Self> {
        let subtile_width = width.div_euclid(usize::try_from(factor)?);
        Ok(Self {
            subtiles: Vec::new(),
            tile_width: width,
            subtile_width,
            invalid_count: 0,
        })
    }

    /// Make all the subtiles for a DEM tile.
    pub fn make_all_subtiles(&mut self, tile: &srtm_reader::Tile) {
        let mut x_offset = 0;
        let mut y_offset = 0;

        while y_offset < self.tile_width {
            if let Some(subtile) = self.make_subtile(tile, x_offset, y_offset) {
                self.subtiles.push(subtile);
            }

            x_offset += self.subtile_width;
            if x_offset >= self.tile_width {
                x_offset = 0;
                y_offset += self.subtile_width;
            }
        }
    }

    /// Make a single subtile.
    fn make_subtile(
        &mut self,
        tile: &srtm_reader::Tile,
        x_offset: usize,
        y_offset: usize,
    ) -> Option<MaxSubTile> {
        let mut max_height = i16::MIN;
        let mut is_valid = false;
        let mut x_end = 0;
        let mut y_end = 0;

        for sub_y in y_offset..(y_offset + self.subtile_width) {
            y_end = sub_y;
            if sub_y >= self.tile_width {
                break;
            }

            for sub_x in x_offset..(x_offset + self.subtile_width) {
                x_end = sub_x;
                if sub_x >= self.tile_width {
                    break;
                }

                let elevation = Self::get_elevation(tile, sub_x, sub_y);
                if elevation > max_height {
                    is_valid = true;
                    max_height = elevation;
                }
            }
        }

        if !is_valid {
            self.invalid_count += 1;
            tracing::error!(
                "No max value found in: {},{}/{},{}",
                tile.longitude,
                tile.latitude,
                x_offset,
                y_offset
            );
            return None;
        }

        #[expect(
            clippy::as_conversions,
            clippy::cast_precision_loss,
            reason = "The inputs are well known"
        )]
        let centre_xy = (
            f32::midpoint(x_offset as f32, x_end as f32),
            f32::midpoint(y_offset as f32, y_end as f32),
        );
        let centre_lonlat = Self::xy_to_latlon(tile, centre_xy);

        let subtile = MaxSubTile {
            lon: centre_lonlat.0,
            lat: centre_lonlat.1,
            max_height: i32::from(max_height),
        };

        Some(subtile)
    }

    /// Get the elevation at the DEM-relative point coordinates of the main tile.
    fn get_elevation(tile: &srtm_reader::Tile, x: usize, y: usize) -> i16 {
        let index = y * tile.resolution.extent() + x;
        #[expect(clippy::panic, reason = "This would be a serious bug")]
        *tile
            .data
            .get(index)
            .unwrap_or_else(|| panic!("Tried to get elevation for point outside of DM: {x}x{y}"))
    }

    /// Convert a DEM-relative x/y coordinate to lon/lat.
    fn xy_to_latlon(tile: &srtm_reader::Tile, (x, y): (f32, f32)) -> (f32, f32) {
        let constant = ARCSEC_PER_DEG / 3.0;

        let lon = f32::from(tile.longitude) + x / constant;
        let lat = f32::from(tile.latitude) + 1.0 - y / constant;
        (lon, lat)
    }

    /// Save the subtiles to disk.
    pub fn save(&self, path: &str) -> Result<()> {
        let bytes: &[u8] = bytemuck::cast_slice(&self.subtiles);
        let mut file = std::fs::File::create(path)?;
        file.write_all(bytes)?;
        Ok(())
    }

    /// Load subtiles from disk.
    pub fn load(path: &str) -> Result<Vec<MaxSubTile>> {
        let bytes = std::fs::read(path)?;

        let slice: &[MaxSubTile] = bytemuck::try_cast_slice(&bytes).map_err(|error| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Error parsing saved file: {error:?}"),
            )
        })?;

        Ok(slice.to_vec())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const SUBTILE_FACTOR: u32 = 10;
    const WIDTH: usize = srtm_reader::Resolution::SRTM3.extent();
    const TOTAL_TILES: usize = (SUBTILE_FACTOR + 1).pow(2) as usize;

    fn make_hgt() -> srtm_reader::Tile {
        let total = WIDTH.pow(2);
        let mut data = vec![0; total];
        data[WIDTH + 2] = 42;
        data[total - 1] = 69;
        srtm_reader::Tile {
            longitude: -2,
            latitude: 54,
            resolution: srtm_reader::Resolution::SRTM3,
            data,
        }
    }

    #[test]
    fn max_subtile_top_left() {
        let tile = make_hgt();
        let mut subtiler = Subtiler::new(SUBTILE_FACTOR, tile.resolution.extent()).unwrap();
        subtiler.make_all_subtiles(&tile);
        assert_eq!(
            subtiler.subtiles[0],
            MaxSubTile {
                lon: -1.9504167,
                lat: 54.950417,
                max_height: 42,
            }
        );
    }

    #[test]
    fn max_subtile_bottom_right() {
        let tile = make_hgt();
        let mut subtiler = Subtiler::new(SUBTILE_FACTOR, tile.resolution.extent()).unwrap();
        subtiler.make_all_subtiles(&tile);
        assert_eq!(subtiler.subtiles.len(), TOTAL_TILES);
        assert_eq!(
            subtiler.subtiles[TOTAL_TILES - 1],
            MaxSubTile {
                lon: -0.99958336,
                lat: 53.999584,
                max_height: 69
            }
        );
    }

    #[test]
    fn save_and_load() {
        let file = tempfile::NamedTempFile::new().unwrap();
        let path = file.path().to_str().unwrap();
        let tile = make_hgt();
        let mut subtiler = Subtiler::new(SUBTILE_FACTOR, tile.resolution.extent()).unwrap();
        subtiler.make_all_subtiles(&tile);
        subtiler.save(path).unwrap();
        let all_subtiles = Subtiler::load(path).unwrap();
        assert_eq!(all_subtiles.len(), TOTAL_TILES);
        assert_eq!(subtiler.subtiles[TOTAL_TILES - 1].max_height, 69i32);
    }
}
