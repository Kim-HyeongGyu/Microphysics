#!/bin/csh -f

rename finite_difference finite_difference0 finite_difference?.png
convert -delay 1 *.png finite_difference.gif
