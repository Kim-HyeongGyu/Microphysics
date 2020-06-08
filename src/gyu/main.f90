program driver
use          namelist_mod
use           file_io_mod, only: write_data
use  model_initialize_mod, only: model_initialize, &
                                 model_close
use   dynamics_driver_mod, only: dynamic_driver
use    physics_driver_mod, only: physics_driver

    integer :: n, k

    !-- Build model 
    call model_initialize()

    !-- Write initial info
    call write_data(1)

    !-- Main
    time_loop: do n = 2, nt
        if ( dyn_adv_scheme == 0 ) then
            continue    ! No dynamics process
        else
            call dynamic_driver(n)
        end if

        if ( phy_adv_scheme == 0 ) then
            continue    ! No physics process
        else
            call physics_driver(n)
        end if

        call write_data(n)
    end do time_loop

    call model_close()

end program driver
