#!/bin/csh -f

rename PPM PPM0 PPM?.png
convert -delay 1 *.png PPM.gif
