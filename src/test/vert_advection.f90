module advection_mod
use            global_mod !only: error_mesg
contains
    subroutine compute_advection(w_full, C, dt, nz, &
                                 dz, scheme, next_C)
!-- Input
! w_full = vertical velocity at full coordinate(nz)
! C      = n-1 time step
! nt     = size of time for iteration
! dt     = length of time step
! nz     = Number of vertical grid
! dz     = depth of model layers
! scheme = differencing scheme, use one of these values:
!           finite_difference = second-order centered finite difference
!           finite_volume     = finite volume method
!           PPM               = piecewise parabolic method
!                               See Colella and Woodward (1984)
!
!-- Output
! next_C = advected quantity
!
! Note! Here, flux form is used for advection term
! FLUX_FORM      = solves for -d(wr)/dt
! ADVECTIVE_FORM = solves for -w*d(r)/dt
!
!       Here, we use Lorenz configuration
!       See Figure 1 in Holdaway et al., (2012) 
!       https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1002/qj.2016

    implicit none
    integer,              intent(in) :: dt, nz
    character(len=*),     intent(in) :: scheme
    real, dimension(nz),  intent(in) :: C, dz
    real, dimension(nz), intent(out) :: next_C
    integer :: ks, ke, kstart, kend
    real    :: wgt   ! weight dz
    real    :: Cwgt  ! weight variable
    real    :: Cdt, cn, Cst
    real    :: dC_dt, zbottom = 0.
    real, dimension(nz)   :: w_full, slope
    real, dimension(nz+1) :: w_half, C_half, flux
    character(len=20)     :: eqn_form = "FLUX_FORM"

    ! vertical indexing
    ks     =    1; ke   = nz
    kstart = ks+1; kend = ke

    ! Make stagged grid for advection
    do k = ks+1, ke
        wgt = dz(k-1) / ( dz(k-1)+dz(k) )
        w_half(k) = w_full(k-1) + wgt*(w_full(k)-w_full(k-1))
    end do
    w_half(1) = 0.; w_half(nz+1) = 0.     ! w = 0 at bottom/top

    ! Set Boundary Condition
    ! most likely w = 0 at these points
    ! w_half(1) = 0.; w_half(nz+1) = 0.     ! Homogeneous Dirichlet BC
    flux(ks) = 0.; flux(ke+1) = 0.          ! Neumann BC

    select case (scheme)
        ! 1) 2nd-order Finite difference scheme {{{
        case ("finite_difference")
            do k = ks+1, ke
                wgt  = dz(k-1) / ( dz(k-1)+dz(k) )
                Cwgt = C(k-1) + wgt*( C(k)-C(k-1) )
                flux(k) = w_half(k)*Cwgt
            end do !}}}

        ! 2) Finite Volume method (FVM) {{{
        case ("finite_volume") 
            ! slope along the z-axis
            call slope_z(C, dz, slope)
            do k = kstart, kend
                if (w_half(k) >= 0.) then
                    if (k == ks) cycle          ! inflow
                    cn  = dt*w_half(k)/dz(k-1)  ! courant number
                    Cst = C(k-1) + 0.5*slope(k-1)*(1.-cn)
                else
                    if (k == ke+1) cycle        ! inflow
                    cn  = -dt*w_half(k)/dz(k)
                    Cst = C(k) - 0.5*slope(k)*(1.-cn)
                end if

                flux(k) = w_half(k) * Cst
                print*, cn

                if (cn > 1.) call error_mesg("Courant number > 1")
            end do
            ! do n = 1, nt-1
            !     do k = 2, nz
            !         C_half(k) = ( dz(k-1)*C(k-1,n) + dz(k)*C(k,n) ) &
            !                   / ( dz(k-1) + dz(k) )
            !     end do
            !
            !     ! Boundary values
            !     C_half(1)    = C( 1,n) - ( C_half( 2)-C( 1,n) )
            !     C_half(nz+1) = C(nz,n) - ( C_half(nz)-C(nz,n) )
            !
            !     flux = w_half*C_half*dt
            !     do k = 1, nz
            !         C(k,n+1) = C(k,n) + ( flux(k)-flux(k+1) ) / dz(k)
            !     end do
            ! end do !}}}

        case default
            call error_mesg("Not setup diff_method option. &
                             please check input.nml")
    end select
    
    ! vertical advective tendency
    select case (eqn_form)
        case ("FLUX_FORM")
            do k = ks, ke
                dC_dt     = - (flux(k+1) - flux(k)) / dz(k)
                next_C(k) = C(k) + dC_dt * dt
            end do
            print*, dC_dt
            ! TODO: SEGMENTATION FAULT
            stop
        case ("ADVECTIVE_FORM")
            do k = ks, ke
                dC_dt     = - ( flux(k+1) - flux(k) )          / dz(k) &
                            - ( C(k)*(w_half(k+1)-w_half(k)) ) / dz(k)
                next_C(k) = C(k) + dC_dt * dt
            end do
        case default
            call error_mesg("No setup equation form.")
    end select
    end subroutine compute_advection


    subroutine slope_z(C, dz, slope, limit, linear)
    real, dimension(nz),  intent(in) :: C, dz
    real, dimension(nz), intent(out) :: slope
    logical,   optional,  intent(in) :: limit, linear

    real    :: grad(2:nz)
    real    :: Cmin, Cmax
    logical :: limiters = .true.
    logical :: dolinear = .true.

    if (present( limit)) limiters = limit
    if (present(linear)) dolinear = linear

    ! compute slope (weighted for unequal levels)
    do k = 2, nz
        grad(k) = (C(k)-C(k-1))/(dz(k)+dz(k-1))
    enddo
    if (dolinear) then
        do k = 2, nz-1
            slope(k) = (grad(k+1)+grad(k))*dz(k)
        enddo
     else
        do k = 2, nz-1
            slope(k) = ( grad(k+1)*(2.*dz(k-1)+dz(k)) + &
                         grad(k  )*(2.*dz(k+1)+dz(k)) ) * dz(k) &
                     / (   dz(k-1) + dz(k) + dz(k+1)  )
         enddo
     endif
     slope(1 ) = 2.*grad(2 )*dz(1 )
     slope(nz) = 2.*grad(nz)*dz(nz)

    ! apply limiters to slope
    if (limiters) then
        do k = 1, nz
            if (k >= 2 .and. k <= n-1) then
                Cmin = min(C(k-1), C(k), C(k+1))
                Cmax = max(C(k-1), C(k), C(k+1))
                slope(k) = sign(1.,slope(k)) *  &
                            min( abs(slope(k)), &
                                2.*(C(k)-Cmin), &
                                2.*(Cmax-C(k))  )
            else
                slope(k) = 0.  ! always slope=0
            endif
        enddo
    endif

    end subroutine slope_z

end module advection_mod
