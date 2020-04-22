#!/bin/csh -f

set compiler = ifort         # "gfortran", "ifort", "pgf90"
set src_dir  = "../src/main"  # direcory path (from current directory)
set FFLAGS   = "-i4 -r8"

set NetCDF_LIB = ${NETCDF}/lib
set NetCDF_INC = ${NETCDF}/include
# set NetCDF_LIB = /usr/local/netcdf/432_intel15/lib
# set NetCDF_INC = /usr/local/netcdf/432_intel15/include


set modules  = (         global.f90 \
                        read_nc.f90 \
                       write_nc.f90 \
                        file_io.f90 \
                     initialize.f90 \
                vert_coordinate.f90 \
                 vert_advection.f90 \
                   microphysics.f90 )

# =========================== Do not touch below code ========================
set execdir = exec

if ( ! -d $execdir ) mkdir $execdir
cd $execdir
\rm ../run.x *.o
foreach mod ($modules)
    ln -sf ../${src_dir}/$mod   # link source code
    $compiler $FFLAGS -c $mod -I $NetCDF_INC   # compile source code
    # $compiler -c $mod -module ./$execdir/ -o ./$execdir/$mod:r.o

    # If you use NETCDF library, add below option
    # -I${NETCDF}/include -L${NETCDF}/lib -lnetcdf
end
ln -sf ../${src_dir}/driver.f90
ln -sf ../input.nml

$compiler $FFLAGS driver.f90 -o ../run.x *.o -L $NetCDF_LIB -lnetcdf -lnetcdff
cd .. 
# ============================================================================

./run.x
