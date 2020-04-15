#!/bin/csh -f

set compiler = ifort         # "gfortran", "ifort", "pgf90"
set FFLAGS   = "-i4 -r8"

set mod = microphysics.f90

# =========================== Do not touch below code ========================
$compiler $FFLAGS -c $mod   # compile source code
$compiler $FFLAGS -o ./run.x dist.f90 *.o
# ============================================================================

./run.x
