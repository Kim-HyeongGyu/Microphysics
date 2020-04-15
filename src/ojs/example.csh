#!/bin/csh -rf

set NetCDF_LIB = /usr/local/netcdf/432_pgi133/lib
set NetCDF_INC = /usr/local/netcdf/432_pgi133/include


pgf90 -c read_nc.f90 -I $NetCDF_INC

pgf90 example.f90 -o ex.exe ./read_nc.o -L $NetCDF_LIB -lnetcdf -lnetcdff

ex.exe
\rm ex.exe *.o
