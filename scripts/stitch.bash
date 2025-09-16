function stitch {
	# Mount Everest
	CENTER_LON=86.925000
	CENTER_LAT=27.988119

	# Caerphilly
	CENTER_LON=-3.3132
	CENTER_LAT=51.6210

	# Cardiff Bay
	CENTER_LON=-3.1230
	CENTER_LAT=51.4898

	# Seattle
	# 47.6061° N, 122.3328° W
	CENTER_LON=-122.3328
	CENTER_LAT=47.6061

	# Utah
	CENTER_LON=-111.7334
	CENTER_LAT=40.2257

	# half-size in metres (600 km box -> half = 300 km)
	HALF_METRES=900000
	# HALF_METRES=50000
	# HALF_METRES=25000

	# resolution in metres (pixel size). set to your desired ground sampling distance
	RES=100

	# input DEM files (can be many .hgt or a VRT). adapt to your files
	# e.g. "srtm_tiles/*.hgt" or an existing merged VRT like "srtm_france.vrt"
	INPUTS="$HOME/Downloads/dems/*.hgt"
	# INPUTS="/publicish/dems/M30/*.hgt"
	INPUTS="/publicish/dems/*.hgt"

	# temporary and final outputs
	VRT="../output/merged.vrt"
	OUT_AEQD="everest_1800km_aeqd_${RES}m.bt"
	OUT_AEQD="caerphilly_100km_aeqd_${RES}m.bt"
	OUT_AEQD="cardiff_50km_aeqd_${RES}m.bt"
	OUT_AEQD="seattle_100km_aeqd_${RES}m.bt"
	OUT_AEQD="seattle_1800km_aeqd_${RES}m.bt"
	OUT_AEQD="utah_1800km_aeqd_${RES}m.bt"

	# --- build a VRT mosaic of all input tiles (fast, non-destructive) ---
	gdalbuildvrt -overwrite "${VRT}" ${INPUTS}

	# --- define AEQD projection centered on Everest ---
	AEQD="+proj=aeqd +lat_0=${CENTER_LAT} +lon_0=${CENTER_LON} +units=m +datum=WGS84 +no_defs"

	echo "Target projection: ${AEQD}"
	# --- compute center coordinates in AEQD (in metres) ---
	# gdaltransform expects lon lat input for EPSG:4326
	read CENTER_X CENTER_Y <<<$(printf "%f %f\n" "${CENTER_LON}" "${CENTER_LAT}" |
		gdaltransform -s_srs EPSG:4326 -t_srs "${AEQD}" |
		awk '{printf "%f %f", $1, $2}')

	echo "AEQD center coordinates: ${CENTER_X}, ${CENTER_Y}"

	# --- compute bounding box in AEQD coords (xmin ymin xmax ymax) ---
	xmin=$(awk -v c="${CENTER_X}" -v h="${HALF_METRES}" 'BEGIN{printf "%.0f", c - h}')
	xmax=$(awk -v c="${CENTER_X}" -v h="${HALF_METRES}" 'BEGIN{printf "%.0f", c + h}')
	ymin=$(awk -v c="${CENTER_Y}" -v h="${HALF_METRES}" 'BEGIN{printf "%.0f", c - h}')
	ymax=$(awk -v c="${CENTER_Y}" -v h="${HALF_METRES}" 'BEGIN{printf "%.0f", c + h}')

	echo "AEQD bbox (m): ${xmin} ${ymin} ${xmax} ${ymax}"

	# --- compute bounding box in lat/lon coords (xmin ymin xmax ymax) ---
	bottom_left=$(echo "-$HALF_METRES -$HALF_METRES" | gdaltransform -s_srs "${AEQD}" -t_srs EPSG:4326)
	ll_xmin=$(echo $bottom_left | cut -d" " -f1)
	ll_ymin=$(echo $bottom_left | cut -d" " -f2)
	top_right=$(echo "$HALF_METRES $HALF_METRES" | gdaltransform -s_srs "${AEQD}" -t_srs EPSG:4326)
	ll_xmax=$(echo $top_right | cut -d" " -f1)
	ll_ymax=$(echo $top_right | cut -d" " -f2)
	echo "Lat/lon bbox (degrees): ${ll_xmin} ${ll_ymin} ${ll_xmax} ${ll_ymax}"

	gdalwarp -overwrite \
		-t_srs "${AEQD}" \
		-te ${xmin} ${ymin} ${xmax} ${ymax} \
		-tr ${RES} ${RES} \
		-r bilinear \
		-of BT \
		-co TILED=YES -co COMPRESS=DEFLATE \
		"${VRT}" "../output/${OUT_AEQD}"

	echo "Wrote AEQD-tiled file: ../output/${OUT_AEQD}"

	# gdal_edit -a_ullr ${ll_xmin} ${ll_ymax} ${ll_xmax} ${ll_ymin} "${OUT_AEQD}"
	# gdal_edit -a_ullr ${CENTER_LON} ${CENTER_LAT} ${CENTER_LON} ${CENTER_LAT} "${OUT_AEQD}"
}
