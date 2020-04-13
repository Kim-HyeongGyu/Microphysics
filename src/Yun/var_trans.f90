    program variable_transformation

    implicit none

    integer :: nz=20, nlev=37
!    integer, intent(in) :: nz
!    integer, intent(in) :: nlev
    real, dimension(nlev) :: qv_in=0.02 ! [kg/kg], specific humidity
    real, dimension(nlev) :: Tin        ! [K]
    real, dimension(nlev) :: Pin        ! 
!    real, dimension(nlev), intent(in) :: qv_in ! [kg/kg], specific humidity
!    real, dimension(nlev), intent(in) :: Tin   ! [K]
!    real, dimension(nlev), intent(in) :: Pin   !
    real, dimension(nz) :: qv_out, Tout, Pout
!    real, dimension(nz), intent(out) :: qv_out, Tout, Pout
    real, dimension(nlev) :: Tv     ! [K], virtual Temperature
    real, parameter :: R  = 287.    ! [J/(kg K)]
    real, parameter :: g  = 9.8     ! [m/s^2]
    real :: H       ! [m] sclae height
!    character(len=*) :: vert_var, temp_var

    qv = 20.
    Tv = Tin*(1+(0.61*qv_in*0.001))
    H = (R*Ts)/g

!    if ( vert_var == p ) then
        
!    endif
!    if ( temp_var == theta ) then

!    endif

    endprogram
