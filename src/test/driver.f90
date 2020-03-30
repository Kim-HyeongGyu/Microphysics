program driver
use          constant_mod, only: show_constant
use   vert_coordinate_mod, only: compute_vert_coord
use         advection_mod, only: compute_advection
implicit none

    ! Set variables
    integer :: nz, it, i 
    integer, parameter :: nt = 100  ! TODO! nt = t_final/dt
    real    :: ztop, zbottom = 0.   ! default bottom [km]
    real    :: GRAVITY, PI
    real, dimension(:), allocatable :: z_full, z_half
    real, dimension(:), allocatable :: T, q, w
    ! real, dimension(n+1) :: array_foreward, array_backward
    character(len=20) :: grid_dz, diff_method
    namelist /setup_list/ nz, ztop, grid_dz, diff_method

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

    T = (/ (I+273, I = 20,1,-1) /)  ! lapse rate 1K/km
    w = (/ (I, I = 1,40,2) /)

    print*, "========= Setup variables ========="
    print*, "Num of z-grid   : ", nz
    print*, "Top of model    : ", ztop
    print*, "Grid type       : ", grid_dz
    print*, "Integration method : ", diff_method
    print*, "==================================="
    ! call show_constant()

    ! Calculate dz
    allocate(z_full(nz), z_half(nz+1))
    call compute_vert_coord(ztop, zbottom, nz, grid_dz, z_full, z_half)

    ! time integration
    ! call compute_advection(w, T, nt, diff_method)

    ! output result

    ! deallocate storage
    if (allocated(z_full)) deallocate(z_full)
    if (allocated(z_half)) deallocate(z_half)

end program driver
