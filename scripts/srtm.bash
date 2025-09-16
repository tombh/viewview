URL="https://www.viewfinderpanoramas.org/Coverage%20map%20viewfinderpanoramas_org3.htm"

function _check_srtm_directory {
	if [ -z "$SRTM3_DIRECTORY" ]; then
		echo "You must set the SRTM3_DIRECTORY environment variable."
		exit 0
	fi
}

# Download all the SRTM3 data for the whole planet. Total of zips is ~16GB
function download_all_srtm3 {
	_check_srtm_directory

	ZIPS="$SRTM3_DIRECTORY/zips"
	mkdir -p "$ZIPS"

	links=$(
		curl -s "$URL" |
			grep -Eoi 'href="[^"]+\.zip"' |
			cut -d'"' -f2
	)

	for link in $links; do
		echo "Parsing $link"
		filename=$(echo "$link" | rev | cut -d'/' -f1 | rev)
		path=$DESTINATION/$filename
		if [ -f "$path" ]; then
			echo "Skipping $filename"
		else
			echo "Downloading $filename"
			curl "$link" --output "$path"
			echo "Sleeping for 15 seconds..."
			sleep 15
		fi
	done
}

# Unzip all the DEM data in a directory. SRTM3 totals to ~70GB.
function unzip_srtms {
	set -eo pipefail
	_check_srtm_directory

	_pushd "$SRTM3_DIRECTORY" || exit
	ZIPS="$SRTM3_DIRECTORY/zips"
	for zip in "$ZIPS"/*.zip; do
		unzip -o "$zip"
		for directory in ./*; do
			[ -d "$directory" ] || continue
			[ "$directory" == "./zips" ] && continue

			find "$directory" -name "*.zip" | while read -r subzip; do
				unzip -o "$subzip"
				rm "$subzip"
			done
			find "$directory" -name "*.hgt" | while read -r hgt; do
				mv "$hgt" .
			done
			rmdir "$directory" || true
		done
	done

	_popd || exit
}
