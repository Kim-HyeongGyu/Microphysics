module vert_coordinate_mod
use global_mod
contains
    subroutine compute_vert_coord(ztop, zbottom, nz, grid_dz, z_full, z_half, dz)
    implicit none
    integer,            intent(in)  :: nz                ! Num of z_full
    real,               intent(in)  :: ztop, zbottom     ! [m]
    character(len=*),   intent(in)  :: grid_dz
    real, dimension(:), allocatable, intent(out) :: z_full, z_half
    real, dimension(:), allocatable, intent(out) :: dz
    integer :: k
    real    :: dzr = 1.05

    if (.not. allocated(    dz)) allocate(    dz(nz  ))
    if (.not. allocated(z_full)) allocate(z_full(nz  ))
    if (.not. allocated(z_half)) allocate(z_half(nz+1))

    z_half(1) = zbottom
    select case (grid_dz)
        case ("constant")
            dz = (ztop - zbottom) / nz
            do k = 1, nz
                z_half(k+1) = z_half(k) + dz(k)
            end do
            z_full = z_half(1:nz) + dz/2.
        case ("stretching")
            dz(1) = ztop*(dzr-1)/(dzr**nz-1)
            do k = 2, nz
                dz(k) = dz(k-1)*dzr
            end do
            ! TODO: optimization
            ! z_half(2:nz+1) = z_half(1:nz) + dz(1:nz)
            ! z_full(1:nz)   = z_half(1:nz) + dz(1:nz)/2.
            do k = 2, nz+1
                ! dz(k) = dz(k-1)*dzr   ! size(dz) /= nz+1
                z_half(k) = z_half(k-1) + dz(k-1)
                z_full(k-1) = z_half(k-1) + dz(k-1)/2.
            end do
        case default
            print*, "Not setup grid_dz option. &
                     please check input.nml"
            stop
    end select

    end subroutine compute_vert_coord

    subroutine interpolate_1d(vert_var, &
                              temp_var, &
                                z_full, &
                                 qv_in, &
                               temp_in, &
                                  w_in, &
                               vert_in, &
                                 P_out, &
                                Th_out, &
                                qv_out, &
                                 w_out, &
                                  Psfc  )

!!!!! IN : character :: vert_var, temp_var
!          integer   :: nz, nlev
!          real, dimension(nlev) :: qv_in, temp_in, vert_in, w_in
!          real, dimension(nz)   :: z_full
!          real, parameter       :: Ps    ! [optional]
!!!! OUT : real, dimension(nz)   :: P_out, T_out, Th_out, qv_out, w_out

    implicit none
    character(len=10),      intent(in)  :: vert_var, temp_var     !! NAMELIST...??
    real, dimension(nlev), intent(in)  :: qv_in      ! [kg/kg]
    real, dimension(nlev), intent(in)  :: temp_in    ! [K]
    real, dimension(nlev), intent(in)  :: w_in       ! [m/s]
    real, dimension(nlev), intent(in)  :: vert_in    ! P [hPa] or Z [m]
    real, dimension(nz),   intent(in)  :: z_full     ! [m]
    real, optional,        intent(in)  :: Psfc  ! [hPa] surface pressure
    real, dimension(:),    intent(out), allocatable :: P_out, Th_out, qv_out, w_out
    real, dimension(nlev) :: z_in, P_in, T_in, Th_in
    real, dimension(nlev) :: Tv         ! [K], virtual temperature
    real, dimension(nlev) :: H          ! [m] scale height
    integer :: i, j
    real    :: d1, d2

    allocate(P_out(nz), Th_out(nz), qv_out(nz), w_out(nz))

!!! :: SETTING VARIABLE for test

    ! qv_in = 0.02
    ! temp_in(1) = 298.
    ! do i = 1, nlev-1
    !     temp_in(i+1) = temp_in(i) - 1.
    ! enddo
    ! vert_in(1) = 1000. 
    ! do i = 1, nlev-1
    !     vert_in(i+1) = vert_in(i) - 25.
    ! enddo 

!!! :: (VARIABLES) to (z, P, T, theta, qv, RH)

    if (present(Psfc)) Ps = Psfc

    if (vert_var == 'p' ) then
        if (temp_var == 'theta' ) then  ! input data : P[hPa] & theta[K]
            Th_in = temp_in
            P_in  = vert_in
            T_in  = Th_in*((P_in/Ps)**(R/Cp))
            Tv    = T_in*(1+(0.61*qv_in))
            H     = (R*Tv)/g
            z_in  = -H*(log(P_in/Ps))
        else                            ! input data : P[hPa] & T[K]
            P_in = vert_in
            T_in = temp_in
            Tv   = T_in*(1+(0.61*qv_in))
            H    = (R*Tv)/g
            z_in = -H*(log(P_in/Ps))
            Th_in = T_in*((Ps/P_in)**(R/cp))
        endif
    else
        if (temp_var == 'theta' ) then  ! input data : z[m] & theta[K]
            print*, ":: INPUT DATA VARIABLE ERROR ::"
            print*, ":: 'Pres' is calculated using 'Temp', and"
            print*, ":: 'Temp' is calculated using 'Pres'."
            print*, ":: Please check the input data variable."
        else                            ! input data : z[m] & T[K]
            z_in = vert_in
            T_in = temp_in
            Tv   = T_in*(1+(0.61*qv_in))
            H    = (R*Tv)/g
            P_in = Ps*exp(-(z_in/H))
            Th_in = T_in*((Ps/P_in)**(R/cp))
        endif
    endif

!!! :: Interpolate to fit nz
    ! TODO: Need Extrapolation
    do i = 1, nz
        do j = 1, nlev
            if ( z_in(j) <= z_full(i) .and. z_full(i) <= z_in(j+1) ) then
                d1 = z_full(i) - z_in(j)
                d2 = z_in(j+1) - z_full(i)
                 P_out(i) =  P_in(j)*(d2/(d1+d2)) +  P_in(j+1)*(d1/(d1+d2))
                 w_out(i) =  w_in(j)*(d2/(d1+d2)) +  w_in(j+1)*(d1/(d1+d2))
                qv_out(i) = qv_in(j)*(d2/(d1+d2)) + qv_in(j+1)*(d1/(d1+d2))
                Th_out(i) = Th_in(j)*(d2/(d1+d2)) + Th_in(j+1)*(d1/(d1+d2))
            endif
        enddo
    enddo

    end subroutine interpolate_1d
end module vert_coordinate_mod
