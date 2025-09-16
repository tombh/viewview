function create_viewshed {
	gdal_viewshed \
		-ox 13210.5375624321 \
		-oy -14580.0434395845 \
		-md 100000 \
		-cc 0 \
		-oz 1.65 \
		caerphilly_100km_aeqd_100m.bt \
		cardiff_viewshed.tiff

	gdal_polygonize \
		-f GeoJSON \
		cardiff_viewshed.tiff \
		cardiff_viewshed.json

	ogr2ogr \
		-f GeoJSON \
		-t_srs 'EPSG:4326' \
		-s_srs '+proj=aeqd +lat_0=51.6210 +lon_0=-3.3132 +units=m +datum=WGS84 +no_defs' \
		cardiff_viewshed_latlon.json \
		cardiff_viewshed.json
}
