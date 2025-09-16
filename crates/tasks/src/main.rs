//! sdf
#![cfg_attr(
    test,
    expect(
        clippy::indexing_slicing,
        clippy::as_conversions,
        clippy::unreadable_literal,
        reason = "Tests aren't so strict"
    )
)]

mod max_subtile;

use color_eyre::Result;
use tracing_subscriber::{Layer as _, layer::SubscriberExt as _, util::SubscriberInitExt as _};

fn main() -> Result<()> {
    let file = "max_subtiles.bin";
    setup_logging()?;

    // let saved = max_subtile::Subtiler::load(file)?;
    // tracing::info!("{}", saved.len());
    // std::process::exit(0);

    let mut subtiler = max_subtile::Subtiler::new(10, srtm_reader::Resolution::SRTM3.extent())?;
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

    subtiler.save(file)?;

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
