module constant_mod
implicit none

    real, parameter :: PI = 3.141592
    real, parameter :: R = 287        ! [J kg-1 K-1]
    real, parameter :: gravity = 9.8  ! [m s-2]

contains
    subroutine show_constant()
        print*, "PI = ", PI
        print*, "R  = ", R, " [J kg-1 K-1]"
        print*, "g  = ", gravity, " [m s-2]"
    end subroutine show_constant

end module constant_mod
