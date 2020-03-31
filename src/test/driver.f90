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
    real    :: T_sfc = 293.15       ! [K]
    real, dimension(:), allocatable :: z_full, z_half
    real, dimension(:), allocatable :: T, q, w
    real, dimension(:,:), allocatable :: Tout
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
    w = sin( (/ (I, I = 1,40,2) /) / 10. )

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
    allocate(Tout(nt,nz))
    call compute_advection(w, T, T_sfc, nt, nz, z_full, z_half, &
                           diff_method, Tout)

    ! output result
    open(unit = 30, file = "output.txt")
    ! open(unit = 30, file = "output.dat", form='unformatted', &
    !      status = "unknown", access='direct', recl=4*nz) 
    do i = 1, 100
        write(30,*) Tout(i,:)
    end do
    close(30)

    ! deallocate storage
    if (allocated(z_full)) deallocate(z_full)
    if (allocated(z_half)) deallocate(z_half)
    if (allocated(Tout))   deallocate(Tout)

end program driver
