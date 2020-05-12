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
use      microphysics_mod, only: make_bins, conc_dist, &
                                 conc_growth, compute_mass
implicit none

    call read_namelist()

    call read_data_init(nlev, lev, temp_in, qv_in, w_in)


    ! Calculate dz
    call compute_vert_coord(ztop, zbottom, nz, vertical_grid, &
                            z_full, z_half, dz)

    ! interplate 1d
    call interpolate_1d(vert_var, temp_var, z_full, qv_in,    &
                         temp_in, w_in, lev, Pinit, Thinit,   &
                         qinit, winit)

    ! Comupte dt using CFL conditin
    call compute_dt_from_CFL(CFL_condition, dz, winit, nt, dt)
    allocate(Th(nz,nt), q(nz,nt), T(nz,nt), w(nz))
    Th(:,1) = Thinit
    T (:,1) = Th(:,1)*((Pinit(:)/Ps)**(R/Cp))
    q(:,1)  = qinit
    w       = winit
!    q      = 0
!    q(5,1) = 100.
    w       = 1.

    allocate(mass(nbin,nz))
    call make_bins()
    call conc_dist()

    call show_setup_variables()    

    ! Dynamic: time integration
    allocate(dm_dt(nbin,nz), dmb_dt(nbin+1,nz))
    do n = 1, nt-1
        call compute_advection( w, Th(:,n), dt, nz, dz,    &
                                vertical_advect, Th(:,n+1) )
        call compute_advection( w, q(:,n), dt, nz, dz,     &
                                vertical_advect,  q(:,n+1) )
        T(:,n+1) = Th(:,n+1)*((Pinit(:)/Ps)**(R/Cp))    ! Theta[K] to T[K]
        ! TODO: Some Th values are zero. Maybe extratpolation problem.
        T  = 293.15 ! For test [K]
        do k = 1, nz 
            call conc_growth(T(k,n+1), q(k,n+1), Pinit(k), &
                             dm_dt(:,k), dmb_dt(:,k))
        end do
        ! TODO: Make code for mass
        ! call compute_mass()
    end do

    call write_data()
 
    print*, "Successfully run!"

end program driver
