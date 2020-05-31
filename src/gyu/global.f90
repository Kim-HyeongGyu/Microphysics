module global_mod
! Global module
! all variable units are MKS based.
implicit none
 
    ! Setup variables
    ! From INPUT/{sounding data}
    real, dimension(:), allocatable :: height_in    ! [m] Height / [hPa] [Pa] pressure
    real, dimension(:), allocatable :: temp_in      ! [K] (potential) temperature
    real, dimension(:), allocatable :: qv_in        ! [g kg] [kg kg] mixing ratio
    real, dimension(:), allocatable :: w_in         ! [m s-1] vertical wind

    ! 1D Model Coordinate variables
    real, dimension(:), allocatable :: z_full       ! Full         (size = nz)
    real, dimension(:), allocatable :: z_half       ! Half(staged) (size = nz+1)
    real, dimension(:), allocatable :: dz           ! Grid length

    ! Interpolated variables
    real, dimension(:), allocatable :: Prs          ! [hPa] pressure
    real, dimension(:), allocatable :: T            ! [K] temperature
    real, dimension(:), allocatable :: THETA        ! [K] potential temperature
    real, dimension(:), allocatable :: qv           ! [kg kg] mixing ratio
    real, dimension(:), allocatable :: W            ! [m s-1] vertical wind

    ! Variables for bins
    real, dimension(  :), allocatable :: radius           ! [m] droplet radius
    real, dimension(  :), allocatable :: radius_boundary  ! [m] radius at boundary
    real, dimension(:,:), allocatable :: mass             ! [kg] droplet mass
    real, dimension(:,:), allocatable :: mass_boundary    ! [kg] mass at boundary
    real, dimension(:,:), allocatable :: Nr               ! [# m-3] number of droplet

    real, dimension(:,:), allocatable :: dm_dt      ! [kg s-1] Mass tendency
    real, dimension(:,:), allocatable :: dmb_dt     ! [kg s-1] at boundary

end module global_mod
