require 'pg'

def make_query(sql)
  conn = PG.connect(dbname: 'world')
  conn.exec(sql) do |result|
    result.each do |row|
      puts row
    end
  end
end

def size_of_array
  dem_width_in_km = 1000
  byte_size = 4
  d1 = (dem_width_in_km * 1000) / 30
  d2 = d1 ** 2
  d3 = d1 ** 3
  bos = (d1.to_f * byte_size) / 3 / 1000
  dem = (d2.to_f * byte_size) / (1000 ** 2)
  size3 = (((d2.to_f / 8) * byte_size) * 2 * 20) / (1000 ** 2)
  size4 = (((d3.to_f / 8) * byte_size) / (1000 ** 3)) * 180
  seconds = (size4.to_f + size3.to_f) / 300
  hours = seconds / 360
  puts "#{bos}kb per band of sight"
  puts "#{dem}MB per DEM array"
  puts "#{dem * 8}MB approx working host RAM"
  puts "#{size3}MB for per-sector reserved ring space on device"
  puts "#{size4 / 180}GB for per-sector bands of sight on device"
  puts "#{(size4 / 180) + (size3 / 1000) + (dem * 4 / 1000)}GB working RAM on device"
  puts "#{size4}GB for all bands of sight"
  puts "#{seconds}s (#{hours} hours) transfer for 3d data"
end

#size_of_array; exit

# portishead-benchmark.bt 168x168
# lat: 51.481, lon: -2.769, width: 5000m
# Intel(R) Celeron(R) CPU  N3060  @ 1.60GHz, 1903 MB
# precompute: 3.089s
# compute: 12.700s

#lon = -2.5
#lat = 51.5
# jogja -7.5, 110.5
lat = -7.5
lon = 110.5
#lat = 51.481
#lon = -2.769
# suspension bridge: 51.4549, -2.6278
lon_lat = "#{lon}, #{lat}"
dem_size = 200000
radius = dem_size / 2
point = "ST_MakePoint(#{lon_lat})::geography"
circle = "ST_Buffer(#{point}, #{radius})::geometry"
square = "Box2D(#{circle})"
projected_square = "ST_SetSRID(#{square}, 4326)"

clip = "ST_Clip(rast, #{projected_square})"

dump_path = "/here/dem.gtiff.tmp"
cleaned_path = "/here/dem.gtiff"
projected_path = "/here/dem.bt"
projection_path = "/here/dem.prj"

`rm -f /here/dem.*`

sql = "
SET postgis.gdal_enabled_drivers = 'ENABLE_ALL';
COPY (
  SELECT ST_AsGDALRaster(
    ST_Union(#{clip}),
    'GTIFF'
  )
  FROM elevation
)
TO '#{dump_path}'
WITH BINARY;
"

make_query sql

# Remove the postgres header
`dd bs=25 skip=1 if=#{dump_path} of=#{cleaned_path} 2>&1 > /dev/null`

tile_radius = '-tr 30 30 -tap'
source_projection = '-s_srs EPSG:4326'
custom_projection = "'+proj=aeqd +lat_0=#{lat} +lon_0=#{lon} +elipse=spheroid +units=m +no_defs'"
destination_projection = "-t_srs #{custom_projection}"
projections = "#{source_projection} #{destination_projection}"
sampling = "-r lanczos"
threads = "-multi"
format = "-of BT"

puts `gdalwarp #{threads} #{sampling} #{tile_radius} #{projections} #{format} #{cleaned_path} #{projected_path}`

json = JSON.parse `gdalinfo -json #{projected_path}`

width = json['size'][0]
height = json['size'][1]

if width != height
  min = [width, height].min
  puts "Cropping from #{width}x#{height} to #{min}x#{min}"
  `mv #{projected_path} #{projected_path}.unclipped`
  puts `gdalwarp #{threads} -s_srs #{custom_projection} #{destination_projection} -ts #{min} #{min} #{format} #{projected_path}.unclipped #{projected_path}`
  width = min
  height = min
end

exit
dir = '/home/tombh/Software/tvs'
`cp -f #{projected_path} #{dir}/etc/input/dem.bt`
`cp -f #{projection_path} #{dir}/etc/input/dem.prj`
`cp -f #{projection_path} #{dir}/etc/output/tvs.prj`
cmd = "cd #{dir} && time etc/bin/tvs -dem_width=#{width} -dem_height=#{height}"
puts "precomputing..."
puts `#{cmd} -is_precompute 2>&1`
puts "computing..."
puts `#{cmd} 2>&1`

`gdal_translate -stats ~/Software/tvs/etc/output/tvs.bt ~/Software/tvs/etc/output/tvs.gtiff`

puts "shading..."
`gdaldem color-relief ~/Software/tvs/etc/output/tvs.gtiff ~/Downloads/colors.txt ~/Downloads/tvs.gtiff`

exit
puts "creating tiles..."
`rm -rf ~/Downloads/NEWTWS`
`gdal2tiles.py -w none -r bilinear -z 1-11 ~/Downloads/tvs.gtiff ~/Downloads/NEWTWS`
puts "merging..."
`python2 ~/Workspace/viewview/merge_tiletrees.py ~/Downloads/NEWTWS ~/Downloads/N51W002`

puts "done"

