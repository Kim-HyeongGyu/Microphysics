#!/bin/csh -f

#rename Nr_ Nr_000 Nr_?.png
#rename Nr_ Nr_00 Nr_??.png
#rename Nr_ Nr_0 Nr_???.png
convert -delay 0.5 *.png Nr.gif
