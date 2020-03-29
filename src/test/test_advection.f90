module advection_mod
contains
    subroutine compute_advection(w, C, diff_method, advection)
    ! Vertical advection = w * (dC/dz)
    ! diff_method: 1. Foward/Centerd finite difference
    !              2. Finite volume method
    implicit none
    real, intent(in)    :: w, C, dt
    real, intent(out), dimension(nz) :: advection

    select case (diff_method)
        case ("finite_difference")
            !  C(n+1,k) - C(n,k)           C(n,k+1) - C(n,k-1)
            ! ------------------- = -w(k) ---------------------
            !          dt                  z(n,k+1) - z(n,k-1)
            C(n+1,k) = C(n,k) - w(n,k)*dt*( C(n,k+1) - C(n,k-1) ) &
                                         /( z(n,k+1) - z(n,k-1) )
        case ("finite_volume")
            print*, "TODO! need to work..."
            stop
        case default
            print*, "Not setup diff_method option. &
                     please check input.nml"
    end select

    end subroutine compute_vert_coord
end module advection_mod
