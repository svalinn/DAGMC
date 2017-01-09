#!/bin/bash

# scans the current directory and for each new file not already existing the
# gallery.rst file, creates a new entry, the person adding the file will however
# need to add descriptive text

file=gallery.rst

# function to write a gallery entry
function gallery_entry()
{
    echo "..  image::" $1 >> $2
    echo "    :alt:   Dummy Text Lorem Ipsum" >> $2
    echo "" >> $2
}

# for each png file
for i in *.png ; do
    if ! grep -q $i gallery.rst ; then
       gallery_entry $i $file
    fi
done
