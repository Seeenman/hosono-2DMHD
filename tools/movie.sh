#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Please specify the name of the input images"
    echo "as well as the name of the output file."
    echo "e.g., movie images/my_image my_output"
    exit
fi


# Set variables
input_pattern="$1%04d.png"  # look for files with endings like 0012.png 0012.png
output_file="./movies/$2.mp4"
frame_rate=10

# Run ffmpeg
ffmpeg -framerate $frame_rate -i $input_pattern -c:v libx264 -pix_fmt yuv420p $output_file
