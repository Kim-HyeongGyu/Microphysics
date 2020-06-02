module substeping_mod

interface time_substeping
    module procedure time_substeping_1d!, time_substeping_2d
end interface

contains
    recursive subroutine time_substeping_1d( W, dz, dt, num_substep )
        implicit none
        real, dimension(:), intent(in)    :: W
        real, dimension(:), intent(in)    :: dz
        real,               intent(inout) :: dt
        integer,            intent(inout) :: num_substep

        real :: courant_number

        if ( size(dz) /= size(W) ) then
            ! Dynamics substeping
            ! Note: size(dz) /= size(W), So we use max value of W for calculate
            ! curant number. Strictly, W must be interpolated as z_full coordinate.
            courant_number = maxval( dt/dz * maxval(abs(W)) )
        else
            ! Physics substeping
            ! courant_number = maxval( dt/dm * abs(dm_dt) )
            courant_number = maxval( dt/dz * abs(W) )
        end if

        if (courant_number > 1) then 
            dt = dt / 2.
            num_substep = num_substep * 2
            call time_substeping_1d( W, dz, dt, num_substep )
        end if
        
    end subroutine time_substeping_1d

end module substeping_mod
