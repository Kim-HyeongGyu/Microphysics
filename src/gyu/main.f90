program driver
use          namelist_mod
use  model_initialize_mod, only: model_initialize
use   dynamics_driver_mod, only: dynamic_driver
use    physics_driver_mod, only: physics_driver

    integer :: n, k

    ! Setup model 
    call model_initialize()

    time_loop: do n = 1, nt-1
        call dynamic_driver()

        call physics_driver()

        ! call write_data()
    end do time_loop

    print*, "Successfully run!"
    ! call model_close()

end program driver
