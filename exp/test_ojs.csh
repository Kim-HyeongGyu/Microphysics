#!/bin/csh -f

set compiler = pgf90         # "gfortran", "ifort", "pgf90"
set src_dir  = "../src/ojs"  # direcory path (from current directory)
set FFLAGS   = "-i4 -r8"

#set NetCDF_LIB = ${NETCDF}/lib
#set NetCDF_INC = ${NETCDF}/include
set NetCDF_LIB = /usr/local/netcdf/432_pgi133/lib
set NetCDF_INC = /usr/local/netcdf/432_pgi133/include


set modules  = (         global.f90 \
                        read_nc.f90 \
                        file_io.f90 )

# =========================== Do not touch below code ========================
set execdir = exec

if ( ! -d $execdir ) mkdir $execdir
cd $execdir

foreach mod ($modules)
    ln -sf ../${src_dir}/$mod   # link source code
    $compiler $FFLAGS -c $mod -I $NetCDF_INC   # compile source code
    # $compiler -c $mod -module ./$execdir/ -o ./$execdir/$mod:r.o

    # If you use NETCDF library, add below option
    # -I${NETCDF}/include -L${NETCDF}/lib -lnetcdf
end
ln -sf ../${src_dir}/test.f90
ln -sf ../input.nml

$compiler $FFLAGS test.f90 -o ../run.x *.o -L $NetCDF_LIB -lnetcdf -lnetcdff
# ============================================================================

../run.x
\rm ../run.x *.o
