module vert_coordinate_mod
contains
    subroutine compute_vert_coord(ztop, nz, grid_dz, z_full, z_half)
    implicit none
    real,    intent(in) :: ztop     ! [km]
    integer, intent(in) :: nz       ! Num of z_full
    character(len=*), intent(in)       :: grid_dz
    real, intent(out), dimension(nz)   :: z_full
    real, intent(out), dimension(nz+1) :: z_half
    integer :: k
    real    :: dz       ! delta z for z_half
    real    :: dzr      ! common ratio for z_half, z_full
    real    :: zbottom  ! [km]

    select case (grid_dz)
        case ("constant_grid")
            zbottom = 0.
            dz = (ztop - zbottom) * 1000. / nz
            z_half(1) = zbottom * 1000.
            do k = 1, nz
                z_half(k+1) = z_half(k) + dz
            end do
            z_full = z_half(1:nz) + dz/2.
        case ("stretching_grid")
            zbottom = 1E-3
            dzr = (ztop*1000.)**(1.0/REAL(nz))
            z_half(1) = zbottom*1000.
            do k = 1, nz
                z_half(k+1) = z_half(k) * dzr
                z_full(k)   = z_half(k) + ((z_half(k+1)-z_half(k))/2.)
            end do
        case default
            print*, "Not setup grid_dz option. &
                     please check input.nml"
            stop
    end select

    end subroutine compute_vert_coord
end module vert_coordinate_mod
