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
│   │   ├ driver.f90           # main code
│   │   ├ global.f90           # variable setting
│   │   ├ file_io.f90          # file I/O
│   │   ├ read_nc.f90          # ├ 
│   │   ├ write_nc.f90         # └ 
│   │   ├ initialize.f90       # variable initialization
│   │   ├ vert_coordinate.f90  # realize vertical coordinate
│   │   ├ vert_advection.f90   # calculate vertical advection
│   │   └ microphysics.f90     # concectration distribution
│   └── user                   # Test directory for User
└── exp
    ├── compile.csh            # compile fortran codes
    ├── input.nml              # namelist
    ├── input                  # For input data (T, w, q, ...)
    ├── output                 # For output data (T, w, q, ...)
    ├── run.x                  # Excute file
    └── exec                   # .mod, .o files

```

## TODO
- [x] Make input/output file formate
- [x] Realize bin structure using 
- [x] Test advection code


