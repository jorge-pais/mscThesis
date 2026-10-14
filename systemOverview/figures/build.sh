#!/bin/sh
# Compile every TikZ figure to PDF, plus a PNG for the notebook to display.
cd "$(dirname "$0")"
for f in *.tikz
do
  jn=$(basename $f .tikz)
  pdflatex $f
  rm $jn.aux $jn.log
  magick -density 250 $jn.pdf -background white -alpha remove -alpha off $jn.png
done
