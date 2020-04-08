#!/bin/csh -f

set compiler = ifort # "gfortran", "ifort", "pgf90"

$compiler -c constant.f90 vert_coordinate.f90 vert_advection.f90
# $compiler -c constant.f90 -module ./{PATH_TO_DIR} -o ./{OUTPUT}
$compiler -o ./a.out driver.f90 *.o
./a.out

rm -f *.mod
rm -f *.o
