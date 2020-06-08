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

        real, dimension(size(W)) :: dz_1
        real :: courant_number

        ! Dynamics & Physics substeping
        ! Note: size(dz) /= size(W), So we use 1st grid at surface
        ! Original grid  : | dz(1) | dz(2) | dz(3) | ... | dz(nz)   |
        ! For substeping : | dz(1) | dz(1) | dz(2) | ... | dz(nz-1) | dz(nz) |
        ! - Dynamic: courant_number = maxval( dt/dz * abs(dz_dt) )
        ! - Physics: courant_number = maxval( dt/dm * abs(dm_dt) )
        dz_1(1)  = dz(1)
        dz_1(2:) = dz(1:)
        courant_number = maxval( dt/dz * abs(W) )

        if (courant_number > 1) then 
            dt = dt / 2.
            num_substep = num_substep * 2
            call time_substeping_1d( W, dz, dt, num_substep )
        end if
        
    end subroutine time_substeping_1d

end module substeping_mod
