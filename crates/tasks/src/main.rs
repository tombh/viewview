//! sdf
#![expect(
    clippy::panic_in_result_fn,
    reason = "This is just code for short tasks, so panicking is better"
)]
#![cfg_attr(
    test,
    expect(
        clippy::indexing_slicing,
        clippy::as_conversions,
        clippy::unreadable_literal,
        reason = "Tests aren't so strict"
    )
)]

mod config;
mod max_subtile;
mod packer;
mod projector;
mod tile;

use clap::Parser as _;
use color_eyre::Result;
use tracing_subscriber::{Layer as _, layer::SubscriberExt as _, util::SubscriberInitExt as _};

use crate::projector::LatLonCoord;

/// The file name for the max subtiles data.
const SUBTILES_FILE: &str = "max_subtiles.bin";

fn main() -> Result<()> {
    setup_logging()?;

    let config = crate::config::Config::parse();
    tracing::info!("Initialising with config: {config:?}",);

    match &config.command {
        crate::config::Commands::Packer(packer_config) => {
            let mut packer = packer::Packer::new(packer_config.clone())?;
            match packer_config.one {
                Some(coordinate) => packer.run_one(LatLonCoord(geo::coord! {
                    x: coordinate.0,
                    y: coordinate.1
                }))?,
                None => packer.run_all()?,
            }
        }
        crate::config::Commands::MaxSubTiles(_) => max_subtiles()?,
    }

    Ok(())
}

/// Setup logging.
fn setup_logging() -> Result<()> {
    let filters = tracing_subscriber::EnvFilter::builder()
        .with_default_directive("info".parse()?)
        .from_env_lossy();
    let filter_layer = tracing_subscriber::fmt::layer().with_filter(filters);
    let tracing_setup = tracing_subscriber::registry().with(filter_layer);
    tracing_setup.init();

    Ok(())
}

/// Create the max subtiles acceleration structure.
fn max_subtiles() -> Result<()> {
    let mut subtiler = max_subtile::Subtiler::new(10, srtm_reader::Resolution::SRTM3.extent());
    let mut count = 0u32;
    let dir = std::path::Path::new("/publicish/dems");
    let entries = std::fs::read_dir(dir)?;
    let total = std::fs::read_dir(dir)?.count();
    for entry in entries {
        let path = entry?.path();
        if path.extension().and_then(|string| string.to_str()) == Some("hgt") {
            tracing::debug!("Loading {path:?}");
            let tile = srtm_reader::Tile::from_file(&path).map_err(|error| {
                color_eyre::eyre::eyre!("Couldn't load SRTM tile ({path:?}): {error:?}")
            })?;
            subtiler.make_all_subtiles(&tile);
            count += 1;
            #[expect(
                clippy::as_conversions,
                clippy::cast_precision_loss,
                reason = "Just for logging"
            )]
            let percentage = (count as f32 / total as f32) * 100.0;
            tracing::info!("{}%", percentage);
        }
    }

    if subtiler.invalid_count > 0 {
        tracing::warn!("{} invalid subtiles", subtiler.invalid_count);
    }

    subtiler.save(SUBTILES_FILE)?;

    Ok(())
}
