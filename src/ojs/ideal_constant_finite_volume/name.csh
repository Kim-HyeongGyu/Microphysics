#!/bin/csh -f

rename finite_volume finite_volume0 finite_volume?.png
convert -delay 1 *.png ideal_finite_volume.gif
