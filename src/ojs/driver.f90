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
                                 conc_growth, compute_conc
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
    dt = 0.01; nt = 1000
    Th(:,1) = Thinit
    T (:,1) = Th(:,1)*((Pinit(:)/Ps)**(R/Cp))
    q(:,1)  = qinit
    w       = winit
    ! q      = 0
    ! q(5,1) = 100.
    ! Th(:,1) = 273.
    ! w       = 1.
    !print*, q(:,1)
    allocate(mass(nbin,nz,nt), mass_boundary(nbin+1,nz,nt))
    call make_bins()
    call conc_dist()

    call show_setup_variables()    

    ! Dynamic: time integration
    allocate(dm_dt(nbin,nz), dmb_dt(nbin+1,nz))
    allocate(drop_num(nbin,nz,nt))
    drop_num = 0
    do k = 1, nz
        drop_num(:,k,1) = Nr
    end do
  

    open(80, file="rb.txt",status="unknown")
    write(80,*) radius_boundary
    open(90,file="Nr.txt",status="unknown")
    do n = 1, nt-1
        call compute_advection( w, Th(:,n), dt, nz, dz,    &
                                vertical_advect,"THETA", Th(:,n+1) )
        call compute_advection( w, q(:,n), dt, nz, dz,     &
                                vertical_advect, "qvapor", q(:,n+1) )
        ! call compute_advection( w, drop_num(:,n), dt, nz, dz,     &
        !                         vertical_advect, "Nc", q(:,n+1) )
        T(:,n+1) = Th(:,n+1)*((Pinit(:)/Ps)**(R/Cp))    ! Theta[K] to T[K]
        
        ! print*, "t=",n
        do k = 1, nz
        ! print*, "z=",k
            call conc_growth(T(k,n+1), q(k,n+1), Pinit(k), &
                             dm_dt(:,k), dmb_dt(:,k))
            ! TODO: test in one layer
            mass(:,k,n+1) = mass(:,k,n) + dm_dt(:,1)*dt
            ! TODO: dqv, dT 
        
        !    dqv(k) = -sum(dm_dt(:,1)*dt)
        !    q(k,n+1)=q(:,n+1)+dqv(k)
        !    dT(k) = -(L*dqv(k))/(rho*Cp)    ! define use L, Cp??
        !    T(k,n+1)=T(:,n+1)+dT(k)

        !print*, dm_dt(:,1)
            call compute_conc(dmb_dt(:,k), drop_num(:,k,n), drop_num(:,k,n+1), &
                              mass(:,k,n),mass(:,k,n+1))
        end do
            do i = 1, nbin
            write(90,*) drop_num(i,1,n)
            end do
    end do
    ! print*, dmb_dt(:,10)


    stop

    call write_data()
 
    print*, "Successfully run!"

end program driver
