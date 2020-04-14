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
            do k = 2, nz+1
                dz(k) = dz(k-1)*dzr
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
                               vert_in, &
                                 T_out, &
                                qv_out, &
                                  Psfc  )

!!!!! IN : character :: vert_var, temp_var
!          integer   :: nz, nlev
!          real, dimension(nlev) :: qv_in, temp_in, vert_in
!          real, dimension(nz)   :: z_full
!          real, parameter       :: Ps    ! [optional]
!!!! OUT : real, dimension(nz)   :: T_out, qv_out

    implicit none
    character(len=*),      intent(in)  :: vert_var, temp_var     !! NAMELIST...??
    real, dimension(nlev), intent(in)  :: qv_in      ! [kg/kg]
    real, dimension(nlev), intent(in)  :: temp_in    ! [K]
    real, dimension(nlev), intent(in)  :: vert_in    ! P [hPa] or Z [m]
    real, dimension(:),    intent(in)  :: z_full     ! [m]
    real, optional,        intent(in)  :: Psfc  ! [hPa] surface pressure
    real, dimension(nz),   intent(out) :: T_out, qv_out
    real, dimension(nz)   :: z_out      ! for test
    real, dimension(nlev) :: z_in, T_in
    real, dimension(nlev) :: Tv         ! [K], virtual temperature
    real, dimension(nlev) :: H          ! [m] scale height
    integer :: i, j
    real    :: d1, d2
    real    :: Ps   ! [hPa] surface pressure

!!!!! :: SETTING VARIABLE

    ! qv_in = 0.02
    ! temp_in(1) = 298.
    ! do i = 1, nlev-1
    !     temp_in(i+1) = temp_in(i) - 1.
    ! enddo
    ! vert_in(1) = 1000. 
    ! do i = 1, nlev-1
    !     vert_in(i+1) = vert_in(i) - 25.
    ! enddo 

!!!!! :: (VARIABLES) to (z, T, qv)
    Ps = 1013.; if (present(Psfc)) Ps = Psfc
    print*, R, g

    if ( vert_var == 'p' ) then
        Tv = temp_in*(1+(0.61*qv_in))
        H = (R*Tv)/g
        z_in = -H*(log(vert_in/Ps)) 
    else
        z_in = vert_in
    endif
    
    if ( temp_var == 'theta' ) then
        if ( vert_var == 'p' ) then 
            T_in = temp_in*((vert_in/Ps)**(R/Cp))
        else
            print*, " :: Without air pressure information,"
            print*, " ::  'Temp' can't be calculated from 'theta'."
            print*, " :: Please Check the 'vertical variable'."
        endif
    else
        T_in = temp_in    
    endif

!!!!! :: Interpolate to fit nz
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

    end subroutine interpolate_1d
end module vert_coordinate_mod
