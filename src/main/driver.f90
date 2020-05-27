program driver
use            global_mod
use           read_nc_mod
use          write_nc_mod
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
    ! 30 min integration -> dt*nt = 30 [min] x 60 [s]
    ! dt = 0.005; nt = 360000
    dt = 0.01 ; nt = 180000
    ! dt = 0.05 ; nt = 36000
    allocate(Th(nz,nt), q(nz,nt), T(nz,nt), w(nz))
    Th(:,1) = Thinit
    T (:,1) = Th(:,1)*((Pinit(:)/Ps)**(R/Cp))
    q (:,1) = qinit
    ! w       = winit
    w       = 0.5

    allocate(mass(nbin,nz,nt), mass_boundary(nbin+1,nz,nt))
    call make_bins()
    call conc_dist()

    call show_setup_variables()    

    ! Dynamic: time integration
    allocate(dm_dt(nbin,nz), dmb_dt(nbin+1,nz))
    allocate(drop_num(nbin,nz,nt))
    drop_num = 0
    ! 1st layer boundary condition
    ! do k = 1, nz
    !     drop_num(:,k,1) = Nr
    ! end do

    allocate(dTemp(nz), dqv(nz))
    do n = 1, nt-1
        if (n*dt == 1200) w=0  ! at 20 min, w = 0
        drop_num(:,1,n) = Nr

        call compute_advection( w, Th(:,n), dt, nz, dz,    &
                                vertical_advect, "THETA", Th(:,n+1), Th(1,1) )
        call compute_advection( w, q(:,n), dt, nz, dz,     &
                                vertical_advect, "qvapor", q(:,n+1), q(1,1) )
        do i = 1, nbin
            call compute_advection( w, drop_num(i,:,n), dt, nz, dz,             &
                                    vertical_advect, "Nc", drop_num(i,:,n+1))!,   &
                                    !drop_num(i,1,1) )
        end do
        T(:,n+1) = Th(:,n+1)*((Pinit(:)/Ps)**(R/Cp))    ! Theta[K] to T[K]

        !do k = 1, nz
        do k = 1, 2
            call conc_growth( T(k,n+1), q(k,n+1), Pinit(k), &
                              dm_dt(:,k), dmb_dt(:,k) )
            mass(:,k,n+1) = mass(:,k,n) + dm_dt(:,k)*dt
            call compute_conc( dmb_dt(:,k), drop_num(:,k,n+1), drop_num(:,k,n+1), &
                               mass(:,k,n), mass(:,k,n+1) )
            ! Online Coupling with T and qv
            dqv(k)    = -sum(dm_dt(:,k)*dt)
            q (k,n+1) =  q(k,n+1)+dqv(k)
            dTemp(k)  = -(L*dqv(k))/(rho*Cp)
            T (k,n+1) =  T(k,n+1)+dTemp(k)
            ! Th(k,n+1) = Th(k,n+1)+dTemp(k)
        end do
        Th(:,n+1) = T(:,n+1)*((Ps/Pinit(:))**(R/Cp))    ! T[K] to Theta[K]
        ! TODO: Check mass...
        qc = sum(drop_num(:,1,n)*mass(:,1,n))
!print*, drop_num(8,:,n)

! print*, (mass(:,1,n)*(3./4.)/pi/rho)**(1./3.)   ! <- radius
! print*, mass(:,1,n+1)
! print*, dmb_dt(:,1)
!print*, n
!print*, drop_num(:,2,n+1)

!    print*, qc, sum(drop_num(:,2,n)), sum(mass(:,2,n))
!if(n==10) stop
!    print*, q(:,n) 
!    print*, drop_num(8,:,n)
print*, qc
if(n==10) stop
    end do
    stop

    call write_data()
 
    print*, "Successfully run!"

end program driver
