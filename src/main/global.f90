module global_mod
implicit none
    
    integer :: n, k, i                ! Iteration
    integer :: num_levels, nz         ! Num of levels
    integer :: nlev
    integer :: nt                     ! time step  [s]
    integer :: dt                     ! delta time [s]
    integer :: t_final                ! Last time  [s]
    integer :: drop_var
    integer :: nbin

    real    :: top_of_atmosphere, ztop = 10000.
    real    :: zbottom = 0.
    real    :: CFL_condition = 0.5
    real    :: rmin, rmax, rratio, Nc, qc
    real    :: Ps  = 1013.      ! [hPa] Surface pressure

    real, parameter :: PI  = 3.141592
    real, parameter :: R   = 287.       ! [J kg-1 K-1]
    real, parameter :: Cp  = 1003.5     ! [J kg-1 K-1] specific heat at constant pressure
    real, parameter :: g   = 9.8        ! [m s-2] gravity 
    real, parameter :: rho = 1000.      ! [kg m-3] water density

    real, dimension(:),     allocatable :: z_full, z_half
    real, dimension(:),     allocatable :: w, dz, Pinit, Thinit, qinit, winit
    real, dimension(:),     allocatable :: lev, temp_in, qv_in, w_in  
    real, dimension(:),     allocatable :: radius, radius_boundary
    ! real, dimension(:,:), allocatable :: mass, mass_boundary 
    real, dimension(:,:,:), allocatable :: mass, mass_boundary 
    real, dimension(:,:),   allocatable :: Th, T, q
    real, dimension(:),     allocatable :: Nr, dN_dlnr
    real, dimension(:,:),   allocatable :: dm_dt, dmb_dt
    real, dimension(:,:,:), allocatable :: drop_num
    
    character(len=20) :: vertical_grid, vertical_advect, mass_scheme
    character(len=10) :: vert_var, temp_var
    character(len=20) :: dist_type
    character(len=10) :: status_case
    logical           :: ventilation_effect

public
contains

    subroutine show_setup_variables()
        print*, "========= Setup variables ========="
        print*, "Num of z-grid      : ", nz
        print*, "Num of time        : ", nt
        print*, "Top of model       : ", ztop, " [m]"
        print*, "Grid type          : ", vertical_grid
        print*, "Integration method : ", vertical_advect
        print*, "==================================="
    end subroutine show_setup_variables

    subroutine show_constant()
        print*, "PI = ", PI
        print*, "R  = ", R, " [J kg-1 K-1]"
        print*, "g  = ", g, " [m s-2]"
    end subroutine show_constant

    subroutine error_mesg(message)
        character(len=*), intent(in) :: message
        print*, "ERROR: ", message
        stop
    end subroutine error_mesg

end module global_mod
