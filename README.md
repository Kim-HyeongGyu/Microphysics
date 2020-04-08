# Microphysics
Repository for Microphysics II class in KNU.

## Running model
Just clone repository and run `./compile.csh` script.
```
$ git clone https://github.com/Kim-HyeongGyu/Microphysics.git
$ cd exp
$ ./compile.csh     # First, set compiler and options.
```

## Structure
```
/microphysics
├── README.md
├── doc    # documentation for box model.
├── src    # source code
│   ├── main
│   │   ├ driver.f90                # code for box model.
│   │   └ ...
│   └── user                        # Test directory for User
│       ├── appendix
│       │   ├── FVM_mit.m
│       │   ├── vert_advection.f90  # advection (in FMS)
│       │   └── vert_advection.html
│       ├── constant.f90            # constant variables
│       ├── driver.f90
│       ├── file_io.f90
│       ├── global.f90              # define variables
│       ├── initialize.f90
│       ├── microphysics.f90
│       ├── vert_advection.f90      # advection (Working on it ..)
│       └── vert_coordinate.f90     # base of Lorenz coordinate
└── exp
    ├── compile.csh        # compile fortran code
    ├── input.nml          # namelist
    ├── INPUT              # For input data (T, w, q, ...)
    ├── plot_output.ncl    # For test
    ├── run.x              # Excute file
    └── exec               # .mod, .o files

```

## TODO
- [ ] Make input/output file formate
- [ ] Realize bin structure using 
- [ ] Test advection code


