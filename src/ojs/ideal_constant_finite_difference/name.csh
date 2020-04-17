#!/bin/csh -f

rename finite_difference finite_difference0 finite_difference?.png
convert -delay 1 *.png ideal_finite_difference.gif
