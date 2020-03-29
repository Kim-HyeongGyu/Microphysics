program driver
use          constant_mod, only: show_constant
use   vert_coordinate_mod, only: compute_vert_coord
! use         advection_mod, only: compute_advection
implicit none

    ! Set variables
    integer :: nz, it
    integer, parameter :: nt = 100  ! TODO! nt = t_final/dt
    real    :: ztop, zbottom
    real    :: GRAVITY, PI
    real, dimension(:), allocatable :: dz
    real, dimension(:), allocatable :: T, q, w
    ! real, dimension(n+1) :: array_foreward, array_backward
    character(len=20) :: grid_dz, diff_method
    namelist /setup_list/ nz, ztop, zbottom, grid_dz, diff_method

    zbottom = 0.    ! default bottom [km]

    open  (unit = 8, file = 'input.nml', delim = 'apostrophe')
    read  (unit = 8, nml  = setup_list) 
    close (unit = 8)

    ! initialization
    ! open  (91,file='./INPUT/T.dat',access='direct',recl=size(ttnd)*8)
    ! read  (91,rec=1) T
    ! close (91)
    ! open  (92,file='./INPUT/q.dat',access='direct',recl=size(ttnd)*8)
    ! read  (92,rec=1) q
    ! close (92)
    ! open  (92,file='./INPUT/w.dat',access='direct',recl=size(ttnd)*8)
    ! read  (92,rec=1) w
    ! close (92)

    print*, "========= Setup variables ========="
    print*, "Num of z-grid   : ", nz
    print*, "Top of model    : ", ztop
    print*, "Bottom of model : ", zbottom
    print*, "Grid type       : ", grid_dz
    print*, "Integration method : ", diff_method
    print*, "==================================="
    ! call show_constant()

    ! Calculate dz
    allocate(dz(nz))
    call compute_vert_coord(ztop, zbottom, nz, grid_dz, dz)    ! [m]

    ! time integration
    ! do it = 1, nt
    !     call compute_advection(w(it,:), T(it,:), diff_method)
    ! end do

    ! output result

    ! deallocate storage
    if (allocated(dz)) deallocate(dz)

end program driver
