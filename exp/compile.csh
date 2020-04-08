#!/bin/csh -f

set compiler = ifort         # "gfortran", "ifort", "pgf90"
set src_dir  = "../src/gyu"  # direcory path (from current directory)
set FFLAGS   = "-i4 -r8"

set modules  = (         global.f90 \
                        file_io.f90 \
                     initialize.f90 \
                vert_coordinate.f90 \
                 vert_advection.f90 \
                   microphysics.f90 )

# =========================== Do not touch below code ========================
set execdir = exec

if ( ! -d $execdir ) mkdir $execdir
cd $execdir

foreach mod ($modules)
    ln -sf ../${src_dir}/$mod   # link source code
    $compiler $FFLAGS -c $mod   # compile source code
    # $compiler -c $mod -module ./$execdir/ -o ./$execdir/$mod:r.o

    # If you use NETCDF library, add below option
    # -I${NETCDF}/include -L${NETCDF}/lib -lnetcdf
end
ln -sf ../${src_dir}/driver.f90

$compiler $FFLAGS -o ../run.x driver.f90 *.o
cd .. 
# ============================================================================

./run.x
