module global_mod
implicit none
    
    integer :: n, k, i                ! Iteration
    integer :: num_levels, nz = 10    ! Num of levels
    integer :: nt                     ! time step  [s]
    integer :: dt                     ! delta time [s]
    integer :: t_final                ! Last time  [s]

    real    :: top_of_atmosphere, ztop = 10000.
    real    :: zbottom = 0.
    real    :: CFL_condition = 0.5

    real, parameter :: PI = 3.141592
    real, parameter :: R = 287        ! [J kg-1 K-1]
    real, parameter :: gravity = 9.8  ! [m s-2]

    real, dimension(:),   allocatable :: z_full, z_half
    real, dimension(:),   allocatable :: w, dz, Tinit, qinit
    real, dimension(:,:), allocatable :: T, q
    character(len=20) :: vertical_grid, vertical_advect

public
contains

    subroutine show_setup_variables()
        print*, "========= Setup variables ========="
        print*, "Num of z-grid      : ", nz
        print*, "Top of model       : ", ztop, " [m]"
        print*, "Grid type          : ", vertical_grid
        print*, "Integration method : ", vertical_advect
        print*, "==================================="
    end subroutine show_setup_variables

    subroutine error_mesg(message)
        character(len=*), intent(in) :: message
        print*, "ERROR: ", message
        stop
    end subroutine error_mesg

end module global_mod