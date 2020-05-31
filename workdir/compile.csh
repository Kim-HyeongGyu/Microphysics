#!/bin/csh -f

set compiler = ifort         # "gfortran", "ifort", "pgf90"
set src_dir  = "../src/gyu"  # direcory path (from current directory)
set FFLAGS   = "-i4 -r8"          # option for ifort
# set FFLAGS   = "-freal-4-real-8"    # option for gfortran


set modules  = (         global.f90 \
                       constant.f90 \
                  error_handler.f90 \
                       namelist.f90 \
                        file_io.f90 \
                vert_coordinate.f90 \
                   microphysics.f90 \
               model_initialize.f90 \
                dynamics_driver.f90 \
                 physics_driver.f90 )

# =========================== Do not touch below code ========================
set execdir = exec

if ( ! -d $execdir ) mkdir $execdir
cd $execdir
\rm ../run.x *.o
foreach mod ($modules)
    ln -sf ../${src_dir}/$mod   # link source code
    $compiler $FFLAGS -c $mod   # compile source code
    # $compiler -c $mod -module ./$execdir/ -o ./$execdir/$mod:r.o

end
ln -sf ../${src_dir}/main.f90

$compiler $FFLAGS main.f90 -o ../run.x *.o
cd .. 

./run.x
# ============================================================================
