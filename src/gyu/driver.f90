program driver
use            global_mod
use           file_io_mod, only: read_namelist, &
                                     read_data, &
                                    write_data
use        initialize_mod, only: compute_dt_from_CFL, &
                                      initialize_end
use   vert_coordinate_mod, only: compute_vert_coord
use         advection_mod, only: compute_advection
implicit none

    call read_namelist()
    call show_setup_variables()

    ! Calculate dz
    call compute_vert_coord(ztop, zbottom, nz, vertical_grid, &
                            z_full, z_half, dz)

    call read_data()
    ! print*, T(:,1), q(:,1)

    ! Comupte dt using CFL conditin
    call compute_dt_from_CFL(CFL_condition, dz, w, nt, dt)

    ! Dynamic: time integration
    do n = 1, nt-1
        call compute_advection(w, T(:,n), dt, nz, dz, &
                               vertical_advect, T(:,n+1))
        ! call compute_advection(w, q(:,n), dt, nz, dz, &
        !                        vertical_advect, q(:,n+1))
    end do

    ! Test conservation quantity
    ! do i = 1, 10
    ! ! do i = 1, nt
    !     ! do k = 1, nz 
    !     !     print*, q(k), qout(k,i)
    !     ! end do
    !     print*, sum(q(:)), sum(qout(:,i))
    !     ! if (i == 1) exit
    ! end do

    call write_data()

    call initialize_end()   ! deallocate storage

    print*, "Successfully run!"

end program driver
