How to add new images to the slideshow
======================================

1. Copy the images you wish to add to the gallery to this folder. The images
   must be in `.png` format.
2. Run the `make_gallery.sh` script. This will append the new image data to
   `gallery.rst`. Alternatively, edit `gallery.rst` yourself if you don't want
   the new picture to be at the end of the slideshow.
3. Edit the alt text for the new images you have added in `gallery.rst`; these
   will be at the bottom of the file. Otherwise they will be captioned with
   Lorem Ipsum.
4. Run the `make_html_entry.sh` script. This will rebuild the html part of the
   gallery automatically.
5. Commit the `.png` file and the the modified `gallery.rst` and
   `slideshow.html` files.
