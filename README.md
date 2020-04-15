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
│   │   ├ constant.f90
│   │   ├ driver.f90           # main code
│   │   ├ file_io.f90
│   │   ├ read_nc.f90
│   │   ├ write_nc.f90
│   │   ├ global.f90
│   │   ├ initialize.f90
│   │   ├ microphysics.f90 
│   │   ├ vert_advection.f90
│   │   └ vert_coordinate.f90       
│   └── user                   # Test directory for User
└── exp
    ├── compile.csh            # compile fortran code
    ├── input.nml              # namelist
    ├── input                  # For input data (T, w, q, ...)
    ├── output                 # For input data (T, w, q, ...)
    ├── run.x                  # Excute file
    └── exec                   # .mod, .o files

```

## TODO
- [x] Make input/output file formate
- [x] Realize bin structure using 
- [x] Test advection code


