module vert_coordinate_mod
use           constant_mod
use           namelist_mod
use      error_handler_mod, only: error_mesg
contains

    subroutine compute_vert_coord(ztop, zbottom, nz, grid_dz, z_full, z_half, dz)
    implicit none
    integer, intent(in)  :: nz                ! Num of z_full
    real,    intent(in)  :: ztop, zbottom     ! [m]
    integer, intent(in)  :: grid_dz
    real, dimension(:), allocatable, intent(out) :: z_full, z_half
    real, dimension(:), allocatable, intent(out) :: dz
    integer :: k
    real    :: dzr = 1.05

    if (.not. allocated(    dz)) allocate(    dz(nz  ))
    if (.not. allocated(z_full)) allocate(z_full(nz  ))
    if (.not. allocated(z_half)) allocate(z_half(nz+1))

    z_half(1) = zbottom
    ! Set vertical grid
    select case (grid_dz)
        case (1)    ! Constant grid
            dz = (ztop - zbottom) / nz
            do k = 1, nz
                z_half(k+1) = z_half(k) + dz(k)
            end do
            z_full = z_half(1:nz) + dz/2.
        case (2)    ! Stretching grid
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
            call error_mesg("Not setup grid_dz option. &
                             please check input.nml")
    end select

    end subroutine compute_vert_coord

    subroutine interpolate_1d( vert_in, &
                               temp_in, &
                                 qv_in, &
                                  w_in, &
                                z_full, &
                                z_half, &
                                 P_out, &
                                 T_out, &
                                Th_out, &
                                qv_out, &
                                 w_out  )

!!!!! IN : character :: vert_var, temp_var
!          integer   :: nz, nlev
!          real, dimension(nlev) :: qv_in, temp_in, w_in, vert_in
!          real, dimension(nz)   :: z_full
!!!! OUT : real, dimension(nz)   :: y_out, Th_out, qv_out, w_out

    implicit none
    real, dimension(:), intent(in) :: vert_in    ! P [hPa] or Z [m]
    real, dimension(:), intent(in) :: qv_in      ! [g kg-1] or [kg kg-1]
    real, dimension(:), intent(in) :: temp_in    ! [K]
    real, dimension(:), intent(in) :: w_in       ! [m s-1]
    real, dimension(:), intent(in) :: z_full     ! [m]
    real, dimension(:), intent(in) :: z_half     ! [m]
    real, dimension(:), intent(out), allocatable :: P_out   ! [hPa] pressure
    real, dimension(:), intent(out), allocatable :: T_out   ! [K] temperature
    real, dimension(:), intent(out), allocatable :: Th_out  ! [K] potential temperature
    real, dimension(:), intent(out), allocatable :: qv_out  ! [kg kg-1] mixing ratio
    real, dimension(:), intent(out), allocatable :: w_out   ! [m s1] vertical wind

    integer                      :: i, j, k
    real                         :: d1, d2
    real                         :: Psfc, Tsfc, qvsfc
    real, dimension(size(qv_in)) :: z_in, P_in, T_in, Th_in
    real, dimension(size(qv_in)) :: Tv            ! [K], virtual temperature
    real, dimension(size(qv_in)) :: H             ! [m] scale height

    if (.not. allocated( P_out)) allocate( P_out(size(z_full)))
    if (.not. allocated( T_out)) allocate( T_out(size(z_full)))
    if (.not. allocated(Th_out)) allocate(Th_out(size(z_full)))
    if (.not. allocated(qv_out)) allocate(qv_out(size(z_full)))
    if (.not. allocated( w_out)) allocate( w_out(size(z_half)))

!!! :: (VARIABLES) to (z, P, T, theta, qv, RH)

    if (surface_data) then
        Psfc  = vert_in(1)
        Tsfc  = temp_in(1)
        qvsfc =   qv_in(1)
    end if

    if (vert_var == 'p' ) then
        if (temp_var == 'theta' ) then  ! input data : P[hPa] & theta[K]
            Th_in = temp_in
            T_in  = Th_in*((P_in/P0)**(R/Cp))
            if (surface_data) T_in(1) = Tsfc
        else                            ! input data : P[hPa] & T[K]
            T_in  = temp_in
            Th_in = T_in*((P0/P_in)**(R/cp))
        endif
        P_in  = vert_in
        Tv    = T_in*(1+(0.61*qv_in))
        H     = (R*Tv)/g
        z_in  = -H*(log(P_in/Psfc))
    else
        if (temp_var == 'theta' ) then  ! input data : z[m] & theta[K]
            print*, "Note! We assumed that Scale height(H) is 8 [km]"
            H     = 8000.
            Th_in = temp_in
            z_in  = vert_in
            if (surface_data) z_in(1) = 0
            P_in  = Psfc*exp(-(z_in/H))
            T_in  = Th_in*( (P_in/P0)**(R/Cp) )
            if (surface_data) T_in(1) = Tsfc
        else                            ! input data : z[m] & T[K]
            T_in = temp_in
            z_in = vert_in
            if (surface_data) z_in(1) = 0
            Tv   = T_in*(1+(0.61*qv_in))
            H    = (R*Tv)/g
            P_in = Psfc*exp(-(z_in/H))
            Th_in = T_in*((P0/P_in)**(R/cp))
        endif
    endif

!!! :: Interpolate & Extrapolation to fit Full coordinate
    do i = 1, size(z_full)
        do j = 1, size(vert_in)

            if ( z_full(i) <= z_in(1) ) then
                d1 = z_in(2) - z_in(1)
                d2 = z_in(1) - z_full(i)
                P_out(i)  =  P_in(1)*((d1+d2)/d1) -  P_in(2)*(d2/d1)
                qv_out(i) = qv_in(1)*((d1+d2)/d1) - qv_in(2)*(d2/d1)
                Th_out(i) = Th_in(1)*((d1+d2)/d1) - Th_in(2)*(d2/d1)

            elseif ( z_in(j) <= z_full(i) .and. z_full(i) <= z_in(j+1) ) then
                d1 = z_full(i) - z_in(j)
                d2 = z_in(j+1) - z_full(i)
                P_out(i)  =  P_in(j)*(d2/(d1+d2)) +  P_in(j+1)*(d1/(d1+d2))
                qv_out(i) = qv_in(j)*(d2/(d1+d2)) + qv_in(j+1)*(d1/(d1+d2))
                Th_out(i) = Th_in(j)*(d2/(d1+d2)) + Th_in(j+1)*(d1/(d1+d2))

            elseif ( z_in(size(vert_in)) <= z_full(i) ) then
                d1 = z_in(size(vert_in)) - z_in(size(vert_in)-1)
                d2 = z_full(i) - z_in(size(vert_in))
                P_out(i)  =  - P_in(size(vert_in)-1) * (     d2/d1) &
                             + P_in(size(vert_in)  ) * ((d1+d2)/d1)
                qv_out(i) = - qv_in(size(vert_in)-1) * (     d2/d1) &
                            + qv_in(size(vert_in)  ) * ((d1+d2)/d1)
                Th_out(i) = - Th_in(size(vert_in)-1) * (     d2/d1) &
                            + Th_in(size(vert_in)  ) * ((d1+d2)/d1)

            endif
        enddo
    enddo

!!! :: Interpolate & Extrapolation to fit Half coordinate
    do i = 1, size(z_half)
        do j = 1, size(vert_in)

            if ( z_half(i) <= z_in(1) ) then
                d1 = z_in(2) - z_in(1)
                d2 = z_in(1) - z_half(i)
                w_out(i)  =  w_in(1)*((d1+d2)/d1) -  w_in(2)*(d2/d1)

            elseif ( z_in(j) <= z_half(i) .and. z_half(i) <= z_in(j+1) ) then
                d1 = z_half(i) - z_in(j)
                d2 = z_in(j+1) - z_half(i)
                w_out(i)  =  w_in(j)*(d2/(d1+d2)) +  w_in(j+1)*(d1/(d1+d2))

            elseif ( z_in(size(vert_in)) <= z_half(i) ) then
                d1 = z_in(size(vert_in)) - z_in(size(vert_in)-1)
                d2 = z_full(i) - z_in(size(vert_in))
                w_out(i)  =  - w_in(size(vert_in)-1) * (     d2/d1) &
                             + w_in(size(vert_in)  ) * ((d1+d2)/d1)

            endif
        enddo
    enddo
    T_out  = Th_out*( (P_out/P0)**(R/Cp) )

    end subroutine interpolate_1d
end module vert_coordinate_mod
