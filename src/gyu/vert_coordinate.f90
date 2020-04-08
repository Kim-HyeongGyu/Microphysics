module vert_coordinate_mod
use            global_mod, only: error_mesg
contains
    subroutine compute_vert_coord(ztop, zbottom, nz,     &
                                  grid_dz, z_full, z_half, dz)
    implicit none
    integer,            intent(in)  :: nz                ! Num of z_full
    real,               intent(in)  :: ztop, zbottom     ! [m]
    character(len=*),   intent(in)  :: grid_dz
    real, dimension(:), allocatable, intent(out) :: z_full, z_half
    real, dimension(:), allocatable, intent(out) :: dz
    integer :: k
    real    :: dzr = 1.05

    if (.not. allocated(    dz)) allocate(    dz(nz  ))
    if (.not. allocated(z_full)) allocate(z_full(nz  ))
    if (.not. allocated(z_half)) allocate(z_half(nz+1))

    z_half(1) = zbottom
    select case (grid_dz)
        case ("constant")
            dz = (ztop - zbottom) / nz
            do k = 1, nz
                z_half(k+1) = z_half(k) + dz(k)
            end do
            z_full = z_half(1:nz) + dz/2.
        case ("stretching")
            dz(1) = ztop*(dzr-1)/(dzr**nz-1)
            do k = 2, nz+1
                dz(k) = dz(k-1)*dzr
                z_half(k) = z_half(k-1) + dz(k-1)
                z_full(k-1) = z_half(k-1) + dz(k-1)/2.
            end do
        case default
            call error_mesg("Not setup grid_dz option. &
                             please check input.nml")
    end select

    end subroutine compute_vert_coord
end module vert_coordinate_mod
