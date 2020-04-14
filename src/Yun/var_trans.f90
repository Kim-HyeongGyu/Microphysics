    program variable_transformation

    implicit none

!!!!! IN : character :: vert_var, temp_var
!          integer   :: nz, nlev
!          real, dimension(nlev) :: qv_in, temp_in, vert_in
!          real, dimension(nz)   :: z_full
!          real, parameter       :: Ps    ! [optional]
!!!! OUT : real, dimension(nz)   :: T_out, qv_out

!    character(len=*), intent(in) :: vert_var, temp_var     !! NAMELIST...??
    integer :: nz=20, nlev=40, i, j
!    integer, intent(in) :: nz, nlev
    real, dimension(nlev) :: qv_in      ! [kg/kg], specific humidity
    real, dimension(nlev) :: temp_in    ! [K]
    real, dimension(nlev) :: vert_in
!    real, dimension(nlev), intent(in) :: qv_in, temp_in, vert_in
    real, dimension(nz+1) :: z_half             !! ***
    real, dimension(nz)   :: z_full, dz, z_out  !! ***
!    real, dimension(nz), intent(in) :: z_full
    real, parameter :: Ps = 1013.  ! [hPa] surface pressure
!    real, parameter, intent(in) :: Ps  ! [hPa] surface pressure
    real, dimension(nz) :: T_out, qv_out
!    real, dimension(nz), intent(out) :: T_out, qv_out
    real, dimension(nlev) :: z_in, T_in
    real, dimension(nlev) :: Tv     ! [K], virtual temperature
    real, dimension(nlev) :: H      ! [m] scale height
    real, parameter :: g  = 9.8     ! [m s-2]
    real, parameter :: R  = 287.    ! [J kg-1 K-1]
    real, parameter :: Cp = 1003.5  ! [J kg-1 K-1] specific heat at constant pressure
    real :: d1, d2
    real :: ztop=16000., dzr=1.05           !! ***

!!!!! :: SETTING VARIABLE

    qv_in = 0.02
    temp_in(1) = 298.
    do i = 1, nlev-1
        temp_in(i+1) = temp_in(i) - 1.
    enddo
    vert_in(1) = 1000. 
    do i = 1, nlev-1
        vert_in(i+1) = vert_in(i) - 25.
    enddo 
    
    dz(1) = ztop*(dzr-1)/(dzr**nz-1)
    do i = 2, nz+1
        dz(i) = dz(i-1)*dzr
        z_half(i) = z_half(i-1) + dz(i-1)
        z_full(i-1) = z_half(i-1) + dz(i-1)/2.
    enddo

!!!!! :: (VARIABLES) to (z, T, qv)

!    if ( vert_var == 'p' ) then
        Tv = temp_in*(1+(0.61*qv_in))
        H = (R*Tv)/g
        z_in = -H*(log(vert_in/Ps)) 
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
    print*, ''
    do i = 1, nz
        do j = 1, nlev
            if ( z_in(j) <= z_full(i) .and. z_full(i) <= z_in(j+1) ) then
                d1 = z_full(i) - z_in(j)
                d2 = z_in(j+1) - z_full(i)
                z_out(i)  =  z_in(j)*(d2/(d1+d2)) +  z_in(j+1)*(d1/(d1+d2))
                T_out(i)  =  T_in(j)*(d2/(d1+d2)) +  T_in(j+1)*(d1/(d1+d2))
                qv_out(i) = qv_in(j)*(d2/(d1+d2)) + qv_in(j+1)*(d1/(d1+d2))
            endif
        enddo
        print*, 'i : ', i, 'z_out : ', z_out(i), 'T_out : ', T_out(i)
    enddo
    print*, ''
    print*, '-----------------------------------------------------------------'
    print*, ''
    print*, 'z_full'
    print*, z_full
    print*, ''
    print*, 'z_in' 
    print*, z_in
    print*, ''
    print*, 'z_out : '
    print*, z_out
    print*, ''
    print*, '-----------------------------------------------------------------'
    print*, ''
    print*, 'T_in'
    print*, T_in
    print*, ''
    print*, 'T_out : '
    print*, T_out

    endprogram
