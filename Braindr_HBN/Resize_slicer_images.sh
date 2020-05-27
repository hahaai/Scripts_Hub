#  Re-size the figures from slicer to fit Braindr display. 


# convert png to jpeg
#magick convert rose.jpg rose.png


rm Crop*
for i in $(find Braindr_images_jpg -iname '*_A*' ! -iname '*Crop*');do
	i_new=${i/braindr_images_jpg/braindr_images_upload}
  echo $i
  convert $i      -resize 160x165^ \
          -gravity center -background black -extent 192x208 $i_new
done


rm Crop*
for i in $(find Braindr_images_jpg -iname '*_S*' ! -iname '*Crop*');do
  echo $i
	i_new=${i/braindr_images_jpg/braindr_images_upload}
  convert $i      -resize 190x190^ \
          -gravity center -background black -extent 192x208  $i_new
done