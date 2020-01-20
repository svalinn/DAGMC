#!/bin/bash

# on the basis of the entries in gallery.rst file
# produce the html snippet that allows to make the
# gallery on the front page

# this function writes the html entry for a given image, as determined
# from the gallery.rst entry
function html_slide_entry()
{
    slide_number="$1"
    number_of_slides="$2"
    imagename="$3"
    alttext="$4"
    output_html="$5"

    echo '  <div class="mySlides fade">' >> $5
    echo '    <div class="numbertext">'$1' / '$2'</div>' >> $5
    echo '    <img src="_images/'$3'" style="width:100%">' >> $5
    echo '    <div class="text">'$4'</div>' >> $5
    echo '  </div>' >> $5
    echo '' >> $5
}

# writes the forward and backward buttons
function html_write_buttons
{
    outfile=$1
    echo '  <a class="prev" onclick="plusSlides(-1)">&#10094;</a>' >> $outfile
    echo '  <a class="next" onclick="plusSlides(1)">&#10095;</a>' >> $outfile
}

# writes the dots for the images
function html_write_dots
{
    dot_type="$1"
    num_dots="$2"
    outfile="$3"
    echo '<div style="text-align:center">' >> $outfile
    for (( i = 1 ; i <= $num_dots ; i++ )) ; do
      echo '  <span class="'$dot_type'" onclick="currentSlide('$i')"></span>' >> $outfile
    done
    echo '</div>' >> $outfile
    echo '' >> $outfile
}


# this is the input rst file which determines, what the html will show
gallery_file="gallery.rst"

# the file we are writing out to
out_file="test.html"

# determine the number of images
num_images=`grep -c ' image::' $gallery_file`

# write the start of the html section
echo '<!-- SLIDESHOW CONTENT -->' > $out_file
echo '<!-- START OF SCRIPTED CONTENT -->' >> $out_file
echo '<br>' >> $out_file
echo '<div class="slideshow-container">' >> $out_file

# loop over the number of slides
for (( i = 1 ; i <= $num_images ; i++ )) ; do
    alttext=`grep ":alt:" $gallery_file | sed -n "$i"p | cut -d ' ' -f6-`
    img_file=`grep "image::" $gallery_file | sed -n "$i"p | awk '{print $3}'`
    echo $alttext $img_file
    html_slide_entry $i $num_images $img_file "$alttext" $out_file
done

# write the buttons
html_write_buttons $out_file

echo '</div>' >> $out_file
echo '<br>' >> $out_file

# write the dots
html_write_dots "dot" $num_images $out_file

echo '<!-- END OF SCRIPTED CONTENT -->' >> $out_file

# now we insert the autogen content into the slideshow file
# determine the line on which we insert the auto get content
breakline=`grep -n "<\!-- WE WILL INSERT AUTO GEN CONENT HERE -->" ../slideshow_empty.html | sed -e s'/:/ /'g | awk '{print $1}'`
breaklineplus1=$(($breakline+1))
filelength=`wc -l ../slideshow_empty.html | awk '{print $1}'`

# first print out the file upto the break point
sed -n 1,"$breakline"p ../slideshow_empty.html > ../slideshow.html
# dump the contents of the autogen file into the slideshow
cat $out_file >> ../slideshow.html
# now print out the rest of the slideshow
sed -n "$breaklineplus1","$filelength"p ../slideshow_empty.html >> ../slideshow.html

# all done :)
