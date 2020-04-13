    program variable_transformation

    implicit none

!    character(len=*), intent(in) :: vert_var, temp_var     !! NAMELIST...??
    integer :: nz=20, nlev=40, I
!    integer, intent(in) :: nz, nlev
    real, dimension(nlev) :: qv_in      ! [kg/kg], specific humidity
    real, dimension(nlev) :: temp_in    ! [K]
    real, dimension(nlev) :: vert_in
!    real, dimension(nlev), intent(in) :: qv_in, temp_in, vert_in
    real, parameter :: Ps = 1013.  ! [hPa] surface pressure
!    real, parameter, intent(in) :: Ps  ! [hPa] surface pressure
    real, dimension(nz) :: z_out, Tout, qv_out
!    real, dimension(nz), intent(out) :: z_out, Tout, qv_out
    real, dimension(nlev) :: z_in, T_in
    real, dimension(nlev) :: Tv     ! [K], virtual temperature
    real, dimension(nlev) :: H      ! [m] sclae height
    real, parameter :: g  = 9.8     ! [m s-2]
    real, parameter :: R  = 287.    ! [J kg-1 K-1]
    real, parameter :: Cp = 1003.5  ! [J kg-1 K-1] specific heat at constant pressure

!!!!! :: SETTING VARIABLE

    qv_in = 0.02
    temp_in = 298.
    vert_in(1) = 1000. 
    do i = 1, nlev-1
        vert_in(i+1) = vert_in(i) - 25.
    enddo 

!!!!! :: (VARIABLES) to (z, T, qv)

!    if ( vert_var == 'p' ) then
        Tv = temp_in*(1+(0.61*qv_in))
        H = (R*Tv)/g
        z_in = -H*(log(vert_in/Ps)) 
        print*, z_in
!    else
!       z_in = vert_in
!    endif

!    if ( temp_var == 'theta' ) then
!        if ( vert_var == 'p' ) then 
!            T_in = temp_in*((vert_in/Ps)**(R/Cp))
!        else
!            print*, " :: Without air pressure information,"
!            print*, " ::  'Temp' can't be calculated from 'theta'."
!            print*, " :: Please Check the 'vertical variable'."
!        endif
!    else
       T_in = temp_in    
!    endif

!!!!! :: Interpolate to fit nz


    endprogram
