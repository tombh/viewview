#!/bin/bash

# Download the raw SRTM elevation data from NASA.

# This data covers most of the world apart from the extremeties of north and south (cite).
# Its resolution is 1° or approximately 30m². As degrees diverge from metres the further
# away from the equator you go, the resolution in fact changes, perhaps to around 60m²
# near the UK for example.

# Inspired by code from the Skadi project by Valhalla, under the MIT license:
# https://github.com/valhalla/skadi/blob/master/scripts/get_data.sh

# Note that when this is run directly on mounted storage like Google's buckets, some operations are
# glacially slow. In fact, it is probably best to only run the final `curl` that downloads the
# files on a mounted bucket. Fortunately, once the `srtmgl1.003.*` files are complete, only
# the `curl` section will in fact run.

# You will need to provide a cookie value for USGS_COOKIE to do the actual dowwnloading. This can
# be done by logging into your Earth Data account at: https://urs.earthdata.nasa.gov and then
# navigating to the $USGS_BASE URL to generate the cookie. It seems that the cookie is linked to
# an IP, so you may need to use something like elinks to retrieve the cooke if you are downloading
# on a remote VM.

# The total size of all zipped files is 91.71GB

USGS_BASE='http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL1.003/2000.02.11/'
PARALLEL_DOWNLOADS=10

# Grab the list of the files
if [ ! -e srtmgl1.003.html ]; then
  echo "$(date): fetching srtm file list"
  curl $USGS_BASE -Ls -o srtmgl1.003.html
fi

# Grab the abbreviated list of hgt files
grep -F '.hgt.zip<' srtmgl1.003.html | sed -e 's@.*href="@@g' -e 's/">.*//g' > srtmgl1.003.list

# Filter out files that are already on the disk (so we don't re-download them)
echo -n > srtmgl1.003.urls
for f in $(cat srtmgl1.003.list); do
  # No zip file
  if [ ! -e $f ]; then
    # No hgt file
    hgt=$(echo $f | sed -e 's/\..*/.hgt/g')
    if [ ! -e $hgt ]; then
      echo "${USGS_BASE}${f}" >> srtmgl1.003.urls
    fi
  fi
done

# Get them onto disk
if [ $(wc -l srtmgl1.003.urls | awk '{print $1}') -gt 0 ]; then
  # Download the zip files
  echo "$(date): downloading SRTM files"
  cat srtmgl1.003.urls | xargs -n1 -P$(PARALLEL_DOWNLOADS) curl -H "Cookie: DATA=$USGS_COOKIE" -s -O
fi

