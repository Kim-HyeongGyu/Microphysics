module vert_coordinate_mod
contains
    subroutine compute_vert_coord(ztop, zbottom, nz, grid_dz, z_full, z_half)
    implicit none
    integer, intent(in) :: nz                ! Num of z_full
    real,    intent(in) :: ztop, zbottom     ! [km]
    character(len=*), intent(in)       :: grid_dz
    real, intent(out), dimension(nz)   :: z_full
    real, intent(out), dimension(nz+1) :: z_half
    integer :: k
    real, dimension(nz) :: dz   ! delta z for z_half

    select case (grid_dz)
        case ("constant_grid")
            dz = (ztop - zbottom) * 1000. / nz
            z_half(1) = zbottom * 1000.
            do k = 1, nz
                z_half(k+1) = z_half(k) + dz(k)
            end do
            z_full = z_half(1:nz) + dz/2.
        case ("stretching_grid")
            print*, "TODO! need to work..."
            stop
            ! do k = 1, nz
            !     dz(k) = dz(k-1)*dzr
            !     print*, dz(k)
            ! end do
            ! calc_ztop = dz(1)*(dzr-1)/(dzr-1)
        case default
            print*, "Not setup grid_dz option. &
                     please check input.nml"
            stop
    end select

    end subroutine compute_vert_coord
end module vert_coordinate_mod
