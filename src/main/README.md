## Step to run box model
1. Edit compiler in compile.csh
2. Run `./compile.csh`

## Namelist exmaple
```
&setup_list
    nz      = 20,                  ! num of grid
    ztop    = 10000,               ! Top of atmosphere [m]
    grid_dz = "constant_grid"      ! "constant_grid", "stretching_grid"
    diff_method = "finite_volume"  ! "finite_difference", "finite_volume"
/
```
