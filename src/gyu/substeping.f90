module substeping_mod

interface time_substeping
    module procedure time_substeping_1d, time_substeping_2d
end interface

contains
    recursive subroutine time_substeping_1d( W, dz, dt, num_substep )
        implicit none
        real, dimension(:), intent(in)    :: W
        real, dimension(:), intent(in)    :: dz
        real,               intent(inout) :: dt
        integer,            intent(inout) :: num_substep

        real :: courant_number

        ! Trick: size(dz) /= size(W), So we use max value of W for calculate
        ! curant number. Strictly, W must be interpolated as z_full coordinate.
        courant_number = maxval( dt/dz * maxval(abs(W)) )

        if (courant_number > 1) then 
            dt = dt / 2.
            num_substep = num_substep * 2
            call time_substeping_1d( W, dz, dt, num_substep )
        end if
        
    end subroutine time_substeping_1d

    !TODO: check input argument
    recursive subroutine time_substeping_2d( dmb_dt, mb, dt, num_substep )
        implicit none
        real, dimension(:,:), intent(in)  :: dmb_dt
        real, dimension(:,:), intent(in)  :: mb
        real,               intent(inout) :: dt
        integer,            intent(inout) :: num_substep

        real :: courant_number
        real, dimension(100) :: dm
        
        dm = mb(2:,1) - mb(1:size(mb)-1,1)

        ! Trick: size(dz) /= size(W), So we use max value of W for calculate
        ! curant number. Strictly, W must be interpolated as z_full coordinate.
        courant_number = maxval( dt/dm * maxval(abs(dmb_dt)) )

        if (courant_number > 1) then 
            dt = dt / 2.
            num_substep = num_substep * 2
            call time_substeping_2d( dmb_dt, mb, dt, num_substep )
        end if

    end subroutine time_substeping_2d
    
end module substeping_mod
