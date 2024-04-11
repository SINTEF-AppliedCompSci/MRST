#/bin/sh
# Put me in 'CO2_Atlas_Norwegian_NorthSea_2011"
mkdir "converted"

for f in */hdr.adf
do
    cd `dirname $f`
    gdal_translate -of AAIGrid hdr.adf "../converted/"`dirname "$f"`
    sed -i 's/-3.4028234663852885981e+38/-1e100/g' "../converted/"`dirname "$f"`
    cd ".."
    
done
