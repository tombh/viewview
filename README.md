# A View Of All Views

This repo is for all the supporting code used to find and display the longest line of sight on the planet.

The main viewshed algorithm is another repo https://github.com/tombh/total-viewsheds

The raw elevation data is currently from (October 2025) https://www.viewfinderpanoramas.org/Coverage%20map%20viewfinderpanoramas_org3.htm Other sources of data are available, notably via AWS's https://registry.opendata.aws/terrain-tiles, but as far as I can tell viewfinderpanoramas.org offers a cleaned version, removing noisy data and filling in voids with other sources.

## Packer

![Map of all the longest line of sight tiles in the world](/assets/world_packed.webp)

This map shows a not-terrible packing of the minimum tiles needed to guarantee searching every line of sight on the planet.

To calculate any [viewshed](https://en.wikipedia.org/wiki/Viewshed), and therefore line of sight, you must inevitably provide more data than ends up being used. This is to say that you don't know the limits of what you can see from a given point until you actually calculate it. The only limit you can calculate beforehand is the longest _theoretical_ line of sight based on the highest points within the region you're interested in.

Here are some examples, they are for 2 peaks of the same height and the maximum distance that they could see each other from:

* 8848m  670km
* 6km    553km
* 4km    452km
* 2km    319km
* 1km    226km
* 0.5km  160km
* 1.65m  9km (height of an average human)

Formula: √(2 * Earth Radius * height) * 2

So as long as you provide the raw data within these theoretical limits then you are at least guaranteed to have complete viewsheds. The worst that can happen is that RAM and CPU cycles are wasted on calculating lines that have already terminated.

We could just cover the world with hundreds of 670km x 670km squares and calculate all the lines of sight inside each one. But that's only really necessary in the Himalayas. But then if we start using different size squares (let's call them tiles), then we face the problem of them not packing well. We start to get overlaps which again introduce lots of wasted RAM and CPU cycles because we're re-calculating regions that have already been done.

So can we strike an optimal balance? This is what the `Packer` in this repo tries to do.

### Steps
1. Create a "lower" resolution version of the global elevation data that for every N degrees/minutes/seconds creates subtile that captures the highest point within itself. So it's lower resolution but it hasn't lost any critical data via the conventional side effects of interpolation. These subtiles are like an accelerating data structure for all the lookups we'll be doing to build the tiles.
2. For each subtile on a popable stack:
  1. Create a tile that fits the highest elevation of all the subtiles it covers.
  2. If the tile overlaps with any other tiles move on to the next subtile.
  3. Remove all the subtiles from the stack that the tile covers.
  4. Repeat this process until no more non-overlapping tiles can be found.
3. Repeat step 2 but allow for overlapping tiles.
4. Once all subtiles are covered, run some cleanup, like removing tiles that are already encompassed by larger tiles.

