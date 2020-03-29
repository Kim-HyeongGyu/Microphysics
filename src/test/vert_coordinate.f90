module vert_coordinate_mod
contains
    subroutine compute_vert_coord(ztop, zbottom, nz, grid_dz, dz)
    implicit none
    integer, intent(in) :: nz                ! [m]
    real, intent(in)    :: ztop, zbottom     ! [km]
    real, intent(out), dimension(nz) :: dz
    character(len=*), intent(in)     :: grid_dz
    integer         :: k
    real, parameter :: s = 1.05
    real :: calc_ztop

    select case (grid_dz)
        case ("constant_grid")
            dz = (ztop - zbottom) * 1000. / nz
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
    end select
    ! if (grid_dz == "constant_grid") then
    !     dz = (ztop - zbottom) * 1000. / nz
    ! end if
    !
    ! if (grid_dz == "stretching_grid") then
    !     
    ! end if
    end subroutine compute_vert_coord
end module vert_coordinate_mod
