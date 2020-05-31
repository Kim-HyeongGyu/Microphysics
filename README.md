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
│   │   ├ read_nc.f90
│   │   ├ write_nc.f90
│   │   ├ initialize.f90       # variable initialization
│   │   ├ vert_coordinate.f90  # realize vertical coordinate
│   │   ├ vert_advection.f90   # calculate vertical advection
│   │   └ microphysics.f90     # concectration distribution
│   └── user                   # Test directory for User
├── exp
│   ├── compile.csh            # compile fortran codes
│   ├── input.nml              # namelist
│   ├── input                  # For input data (T, w, q, ...)
│   ├── output                 # For output data (T, w, q, ...)
│   ├── run.x                  # Excute file
│   └── exec                   # .mod, .o files
└── workdir  # Test for optimizated codes   
```

## What's New?
- Update namelist option for convinience
- *gfortran*, *ifort*, *pgf90*(+2008 version) compilers can use.
- Do not need *NetCDF* library

## TODO
- [ ] Testing physics driver 
- [ ] Porting redistribution from ncl code
- [ ] Advection test in streching grid (Conservation test)
- [ ] +Collection efficiency (from Hall (1980) - JAS)
- [ ] Terminal velocity (from Beard (1977) -MWR)
- [ ] SCE (from Bott (1989) -MWR)

