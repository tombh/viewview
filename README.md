# A View Of All Views

How to use the data from https://www.viewfinderpanoramas.org/Coverage%20map%20viewfinderpanoramas_org3.htm fed into https://github.com/tombh/total-viewsheds and display on views.tombh.co.uk

## Design

How to find a tiling of the world's landmass, where the size of each tile must meet a minimum size based on the elevations that that the tile covers?

* There are approximately 200,000 600km² tiles covering all the land on the planet.
* It'd be good to both use larger tiles, say over Everest, and smaller tiles for islands etc.

1. Create a version of the global DEM data where for every N degree minute/second subtiles I find the highest point. Record the elevation, the highest points lat/lon and the centre of the subtile.
2. For each subtile:
  1. Create an actual tile that fits the highest elevation in all the tiles it covers.
  2. The tile should have a minimum overlap outside any subtiles.
  3. The tile should have a minimum overlap with other tiles.
3. Move on to the next subtile. If it's already been covered by a previous increased tile then move on.
4. Score the final tile based on tile overlaps and going outside the subtiles. Hopefully this is all fast enough that we can brute force a low score. How to change the starting conditions?

### Theoretical longest lines of site based on highest points
These values are for 2 peaks of the given height and the maximum distance that they could see each other:
* 8848m 670km
* 8km   638km
* 7km   597km
* 6km   553km
* 5km   504km
* 4km   452km
* 3km   391km
* 2km   319km
* 1km   226km
* 0.5   160km
* 0     9km

Formula: √(2 * 6,371 * height) * 2

## Assumptions
PostGIS's WGS84 (SRID 4326) projection uses the EGM96 revisions. See: https://en.wikipedia.org/wiki/World_Geodetic_System
_
