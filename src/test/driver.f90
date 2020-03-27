program driver
use constant_mod, only: show_constant
! use         advection_mod, only: compute_advection
! use   vert_coordinate_mod, only: compute_vert_coord
implicit none

    ! Set variables
    integer :: nz
    real    :: ztop
    real    :: GRAVITY, PI
    ! real, dimension(n+1) :: array_foreward, array_backward
    character(len=20) :: grid_dz
    namelist /setup_list/ nz, ztop, grid_dz

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
    
    print*, nz, ztop

    ! dz 계산
    ! dz = compute_dz(ztop, nz, grid_dz)    ! [m]

    ! time integration
    ! do it = 1, nt
    !     call compute_advection(w, T)
    ! end do

    ! output result

end program driver
