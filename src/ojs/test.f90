program driver
use            global_mod
use            read_nc_mod
use           file_io_mod, only: read_namelist, &
                                read_data_init
!                                    write_data
!use        initialize_mod, only: compute_dt_from_CFL, &
!                                      initialize_end
!use   vert_coordinate_mod, only: compute_vert_coord
!use         advection_mod, only: compute_advection
! use      microphysics_mod, only: make_bin
implicit none

    call read_namelist()

    call read_data_init(nlev,lev,temp_in,qv_in,w)


    ! Calculate dz
!    call compute_vert_coord(ztop, zbottom, nz, vertical_grid, &
!                            z_full, z_half, dz)



    print*, "Successfully run!"

end program driver
