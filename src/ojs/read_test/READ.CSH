#!/bin/csh -xvf
#
#set NetCDF_LIB = /usr/local/netcdf/363_pgi716/lib
#set NetCDF_INC = /usr/local/netcdf/363_pgi716/include
set NetCDF_LIB = /usr/local/netcdf/432_pgi133/lib
set NetCDF_INC = /usr/local/netcdf/432_pgi133/include
set FNAME = read_q_nc
#
##--[ Compile ]
#
pgf90 -c $FNAME.f90 -I $NetCDF_INC  
pgf90 -o $FNAME.exe $FNAME.o -L $NetCDF_LIB -lnetcdf -lnetcdff 
#
##--[ Execution ]
#
$FNAME.exe
\rm $FNAME.exe $FNAME.o
#
