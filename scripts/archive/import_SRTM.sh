#!/bin/bash

# Import raw SRTM elevation data into PostGiS database.

# Requires GNU `parallel`: `apt-get install parallel`

# Asssuming Debian Jessie(ish), you'll need something like:
# ```
# apt-get install postgis
# su - postgres
# psql
# create database world;
# \connect world
# create extension postgis;
# ```

# Notes about arguments used:
# -C Constraints like SRID, alignment, resolution, etc. Raises errors if new rasters don't conform.
# -x Stops the -C arg from enforcing extent constraints, so that new rasters can be outside the
#    limits of previous rasters.
# -d DROP TABLE
# -a Append raster to existing table.
# -I Create GiST index (this should be done after import as each call locks the table).
# -Y Use COPY instead of INSERT, it's faster but harder to debug.

# TODO Explore the optimal tile size.

# Master file is for Bristol, UK
MASTER_HGT='N51W003.SRTMGL1.hgt.zip'
# Assume that all the SRTM data was downloaded by our `download_SRTM.sh` script
HGT_LIST=$(cat srtmgl1.003.list)

import(){
  echo "Importing $1..."
  HGT_FILE=~/$(basename $1 .SRTMGL1.hgt.zip).hgt
  gunzip -c $1 > $HGT_FILE
  if [[ "$1" == "$MASTER_HGT" ]]; then
    raster2pgsql -C -x -d $HGT_FILE public.elevation > $HGT_FILE.sql
  else
    raster2pgsql -Y -a $HGT_FILE public.elevation > $HGT_FILE.sql
  fi
  cat $HGT_FILE.sql | psql -d world
  rm $HGT_FILE $HGT_FILE.sql
}
export -f import

# The first file can be used to set constraints for all further imports.
import $MASTER_HGT

for file in $HGT_LIST; do
  if [[ "$file" != "$MASTER_HGT" ]]; then
    sem -j4 import $file
  fi
done

sem -wait

# Leave indexing, vacuuming and analyzing right to the end.
# Indexing takes somewhere around an hour.
# CREATE INDEX ON "public"."elevation" USING gist (st_convexhull("rast"));
# Vacuuming takes about 1.5 hours + ...
# VACUUM ANALYZE "public"."elevation";

