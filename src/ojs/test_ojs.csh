#!/bin/csh -rf

set compiler = pgf90
set FFLAGS = "-i4 -r8"


set NetCDF_LIB = /usr/local/netcdf/432_pgi133/lib
set NetCDF_INC = /usr/local/netcdf/432_pgi133/include

set modules  = (         global.f90 \
                        read_nc.f90 \
                       write_nc.f90 \
                        file_io.f90 \
                     initialize.f90 \
                vert_coordinate.f90 \
                 vert_advection.f90 ) 

set execdir = exec
#echo $modules

if (! -d $execdir ) mkdir $execdir
cd $execdir

foreach mod ( $modules )
  ln -sf ../$mod


     $compiler $FFLAGS  -c $mod -I $NetCDF_INC
end

ln -sf ../driver.f90
ln -sf ../input.nml


$compiler $FFLAGS driver.f90 -o ../ex.exe *.o -L $NetCDF_LIB -lnetcdf -lnetcdff

../ex.exe
\rm ../ex.exe *.o 
