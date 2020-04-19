program driver
use            global_mod
use            read_nc_mod
use           write_nc_mod
use           file_io_mod, only: read_namelist, &
                                read_data_init, &
                                    write_data
use        initialize_mod, only: compute_dt_from_CFL, &
                                      initialize_end
use   vert_coordinate_mod, only: compute_vert_coord , &
                                 interpolate_1d
use         advection_mod, only: compute_advection
use      microphysics_mod, only: make_bins, conc_dist
implicit none

    call read_namelist()

    call read_data_init(nlev, lev, temp_in, qv_in, w_in)


    ! Calculate dz
    call compute_vert_coord(ztop, zbottom, nz, vertical_grid, &
                            z_full, z_half, dz)
    call make_bins()
    call conc_dist()

    ! interplate 1d
    call interpolate_1d(vert_var, temp_var, z_full, qv_in,    &
                         temp_in, w_in, lev, Tinit, qinit, winit)

    ! Comupte dt using CFL conditin
    call compute_dt_from_CFL(CFL_condition, dz, winit, nt, dt)
    allocate(T(nz,nt), q(nz,nt), w(nz) )
    T(:,1) = Tinit
    q(:,1) = qinit
    w = winit
    q      = 0
    q(5,1) = 100.
    w = 1.
    print*, q(:,1)
 
    call show_setup_variables()    


    ! Dynamic: time integration
    do n = 1, nt-1
        call compute_advection(w, T(:,n), dt, nz, dz, &
                               vertical_advect, T(:,n+1))
        call compute_advection(w, q(:,n), dt, nz, dz, &
                               vertical_advect, q(:,n+1))
    end do

    do i = 1, nt
        print*, sum (q(:,i))
    end do

    call write_data()
 
    print*, "Successfully run!"

end program driver
