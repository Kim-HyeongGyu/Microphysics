program driver
use            global_mod
use           file_io_mod, only: read_namelist, &
                                read_data_init, &
                                    write_data
use        initialize_mod, only: compute_dt_from_CFL, &
                                      initialize_end
use   vert_coordinate_mod, only: compute_vert_coord
use         advection_mod, only: compute_advection
use      microphysics_mod, only: make_bins, conc_dist
implicit none

    call read_namelist()

    ! Calculate dz
    call compute_vert_coord(ztop, zbottom, nz, vertical_grid, &
                            z_full, z_half, dz)

    call read_data_init()
    call make_bins()
    call conc_dist()

    ! Comupte dt using CFL conditin
    call compute_dt_from_CFL(CFL_condition, dz, w, nt, dt)
    allocate(T(nz,nt), q(nz,nt))
    T(:,1) = Tinit
    q(:,1) = qinit
    nt = 40!; dt = 10
    call show_setup_variables()

    ! Dynamic: time integration
    do n = 1, nt-1
        call compute_advection(w, T(:,n), dt, nz, dz, &
                               vertical_advect, T(:,n+1))
        call compute_advection(w, q(:,n), dt, nz, dz, &
                               vertical_advect, q(:,n+1))
    end do

    ! Test conservation quantity
    ! do i = 1, nt
    !     print*, sum(q(:,1)), sum(q(:,i))
    !     ! if (i == 1) exit
    ! end do
    call write_data()

    call initialize_end()   ! deallocate storage

    print*, "Successfully run!"

end program driver
