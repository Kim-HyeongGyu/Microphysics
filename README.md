# Microphysics
Repository for Microphysics II class in KNU.

## Structure
```
/microphysics
├ README.md
├ doc    # documentation for box model.
└ src    # source code
    ├ main
    │   ├ driver.f90              # code for box model.
    │   └ microphysics.f90        # Main code (e.g. module, ...)
    └ test
        ├ advection_mod.mod
        ├ a.out
        ├ append_from_fms
        │   ├ vert_advection.f90  # advection (in FMS)
        │   └ vert_advection.html
        ├ compile.csh             # compile fortran code
        ├ constant.f90            # constant variables
        ├ driver.f90
        ├ INPUT                   # input data (T, w, q, ...)
        ├ input.nml               # namelist
        ├ microphysics.f90
        ├ plot_output.ncl         # for test
        ├ vert_advection.f90      # advection (Working on it ..)
        └ vert_coordinate.f90     # base of Lorenz coordinate
```

## Next meeting
- [ ] 2020.04.01 14:00~

