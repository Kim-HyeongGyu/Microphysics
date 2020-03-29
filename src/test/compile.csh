#!/bin/csh -f

set compiler = ifort # "gfortran", "ifort"

$compiler -c constant.f90 vert_coordinate.f90
$compiler -o ./a.out driver.f90 *.o
