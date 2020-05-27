module initialize_mod
use            global_mod
contains

    subroutine compute_dt_from_CFL(courant_number, dz, w, nt, dt)
        implicit none
        real,               intent(in) :: courant_number
        real, dimension(:), intent(in) :: w, dz
        integer,           intent(out) :: nt
        real,              intent(out) :: dt

        dt = minval(courant_number*(dz/abs(w)))
        ! TODO: some compiler makes zero division err
        ! where (w /= 0.)
        !     dt_CFL = CFL_condition*(dz/abs(w))
        ! elsewhere
        !     dt_CFL = maxval(dt_CFL)
        ! end where
        ! dt = minval(CFL)
        nt = int(t_final/dt)
    end subroutine compute_dt_from_CFL

    subroutine initialize_end()
        if (allocated(z_full)) deallocate(z_full)
        if (allocated(z_half)) deallocate(z_half)
        if (allocated(    dz)) deallocate(    dz)
        if (allocated(     T)) deallocate(     T)
        if (allocated(     q)) deallocate(     q)
        if (allocated(     w)) deallocate(     w)
    end subroutine initialize_end

end module initialize_mod
