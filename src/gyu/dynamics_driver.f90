module dynamics_driver_mod
use             global_mod
use           namelist_mod
use           constant_mod
use      error_handler_mod, only: error_mesg
contains

    subroutine dynamic_driver()
        integer :: i
        
        call compute_advection( dt, W, dz, THETA )
        call compute_advection( dt, W, dz, qv    )

        ! Advect each size of droplet in bins
        drop_loop: do i = 1, nbin, 1
            call compute_advection( dt, W, dz, Nr(i,:) )
        end do drop_loop

        ! Convert Theta[K] to T[K] for physics process
        T(:) = THETA(:)*((Prs(:)/P0)**(R/Cp))

    end subroutine dynamic_driver

    subroutine compute_advection( dt, w_half, dz, C )  ! {{{
    !-- Input
    ! dt     = length of time step
    ! w_half = vertical velocity at half coordinate(nz+1)
    ! nt     = size of time for iteration
    ! dz     = depth of model layers
    !
    ! from namelist
    ! dyn_adv_scheme = differencing scheme, use one of these values:
    !  1: finite_difference = second-order centered finite difference
    !  2: finite_volume     = finite volume method
    !  3: PPM               = piecewise parabolic method
    !                        See Colella and Woodward (1984)
    !  4: PPM               = PPM but Lin (2003) limeter used
    !
    !-- Output
    ! C = advected quantity
    !
    ! Note! Here, flux form is used for advection term
    ! FLUX_FORM      = solves for -d(wr)/dt
    ! ADVECTIVE_FORM = solves for -w*d(r)/dt
    !
    ! Here, we use Lorenz configuration
    ! See Figure 1 in Holdaway et al., (2012) 
    ! https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1002/qj.2016

    implicit none
    real,                  intent(in   ) :: dt
    real,    dimension(:), intent(in   ) :: w_half
    real,    dimension(:), intent(in   ) :: dz
    real,    dimension(:), intent(inout) :: C

    integer :: k, kk
    integer :: ks, ke, kstart, kend
    real    :: wgt   ! weight dz
    real    :: Cwgt  ! weight variable
    real    :: dC_dt, zbottom = 0.
    real    :: tt, cn, Csum, dzsum, dtw
    real    :: xx, a, b, Cm, C6, Cst, Cdt
    real, dimension(0:3,size(C)) :: zwt
    real, dimension(size(C))     :: slp, C_left, C_right
    real, dimension(size(C))     :: slope
    real, dimension(size(C)+1)   :: flux
    character(len=20)       :: eqn_form = "FLUX_FORM"
    ! character(len=20)       :: eqn_form = "ADVECTIVE_FORM"
    logical :: linear, test_1
    logical :: do_outflow_bnd = .true.

    ! set default values for optional arguments
    if (dyn_adv_scheme == 1) do_outflow_bnd = .false.

    ! vertical indexing
    ks = 1;  ke = size(C)

    ! start and end indexing for finite volume fluxex
    kstart = ks+1; kend = ke
    if (do_outflow_bnd) then
        kstart = ks;  kend = ke+1
    end if

    ! Make stagged grid for advection
    if ( size(w_half) /= size(C)+1 )then
        call error_mesg("vertical dimension of input arrays inconsistent")
    end if

    ! Set Boundary Condition (Homogeneous Dirichlet BC)
    ! most likely w = 0 at these points
    if (do_outflow_bnd) then
        flux(ks)   = 0.
        flux(ke+1) = 0.
    else
        flux(ks)   = w_half(ks)*C(ks)
        flux(ke+1) = w_half(ke+1)*C(ke)
    end if

    select case (dyn_adv_scheme)
        ! 1) 2nd-order Finite difference scheme {{{
        case (1)
            do k = ks+1, ke
                wgt  = dz(k-1) / ( dz(k-1)+dz(k) )
                Cwgt = C(k-1) + wgt*( C(k)-C(k-1) )
                flux(k) = w_half(k)*Cwgt
            end do !}}}

        ! 2) Finite Volume method (FVM) {{{
        case (2) 
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

                if (cn > 1.) call error_mesg("Courant number > 1")
            end do !}}}

        ! 3) Piecewise Parabolic Method, Colella and Woodward (1984) {{{
        case (3:4)
            zwt = 0 ! Note! Avoid for Nan value occur
            call compute_weights(dz, zwt)
            call slope_z(C, dz, slp, linear=.false.)        ! Equation 1.7
            do k = ks+2, ke-1
                C_left(k) = C(k-1) + zwt(1,k)*(C(k)-C(k-1)) &
                                   - zwt(2,k)*slp(k)        &
                                   + zwt(3,k)*slp(k-1)      ! Equation 1.6 
                C_right(k-1) = C_left(k)
                ! Or, we can use Equation 1.9 (Need condition)
                ! C_rihgt(k) = (7./12.)*(a(k)+a(k+1)) - (1./12.)*(a(k+2)+a(k-1))
                ! coming out of this loop, all we need is r_left and r_right
            enddo

            ! boundary values  ! masks ???????
            C_left (ks+1) = C(ks+1) - 0.5*slp(ks+1)
            C_right(ke-1) = C(ke-1) + 0.5*slp(ke-1)

            ! pure upstream advection near boundary
            ! r_left (ks) = r(ks)
            ! r_right(ks) = r(ks)
            ! r_left (ke) = r(ke)
            ! r_right(ke) = r(ke)

            ! make linear assumption near boundary
            ! NOTE: slope is zero at ks and ks therefore
            !       this reduces to upstream advection near boundary
            C_left (ks) = C(ks) - 0.5*slp(ks)
            C_right(ks) = C(ks) + 0.5*slp(ks)
            C_left (ke) = C(ke) - 0.5*slp(ke)
            C_right(ke) = C(ke) + 0.5*slp(ke)

            if (dyn_adv_scheme == 4) then
                ! limiters from Lin (2003), Equation 6 (relaxed constraint)
                do k = ks, ke
                C_left (k) = C(k) - sign( min(abs(slp(k)),       &
                                          abs(C_left (k)-C(k))), &
                                          slp(k) )  ! (B3)
                C_right(k) = C(k) + sign( min(abs(slp(k)),       &
                                          abs(C_right(k)-C(k))), &
                                          slp(k) )  ! (B4)
                enddo
            else! limiters from Colella and Woodward (1984), Equation 1.10
                do k = ks, ke
                    test_1 = (C_right(k)-C(k))*(C(k)-C_left(k)) <= 0.0
                    if (test_1) then        ! 1.10 (1)
                        C_left (k) = C(k)
                        C_right(k) = C(k)
                    endif
                    if (k == ks .or. k == ke) cycle
                    Cm = C_right(k) - C_left(k)
                    a = Cm*(C(k) - 0.5*(C_right(k) + C_left(k)))
                    b = Cm*Cm/6.
                    if (a >  b) C_left (k) = 3.0*C(k) - 2.0*C_right(k)  ! 1.10 (2)
                    if (a < -b) C_right(k) = 3.0*C(k) - 2.0*C_left (k)  ! 1.10 (3)
                enddo
            endif

            ! compute fluxes at interfaces {{{
            tt = 2./3.
            do k = kstart, kend ! ks+1, nz
                if (w_half(k) >= 0.) then ! w = positive
                    if (k == ks) cycle    ! inflow
                        cn = dt*w_half(k)/dz(k-1)   ! Courant number
                        kk = k-1
                    ! extension for Courant numbers > 1
                    if (cn > 1.) then
                        Csum = 0.; dzsum = 0.
                        dtw  = dt*w_half(k)
                        do while (dzsum+dz(kk) < dtw)
                            if (kk == 1) exit
                            dzsum = dzsum + dz(kk)
                            Csum  =  Csum +  C(kk)
                            kk    =    kk -1
                        enddo
                        xx = (dtw-dzsum)/dz(kk)
                    else
                        xx = cn     ! y = u*dt (1.13)
                    endif
                    Cm = C_right(kk) - C_left(kk)
                    C6 = 6.0*(C(kk) - 0.5*(C_right(kk) + C_left(kk)))   ! (1.5)
                    if (kk == ks) C6 = 0.
                    Cst = C_right(kk) - 0.5*xx*(Cm - (1.0 - tt*xx)*C6)  ! (1.12)
                    ! extension for Courant numbers > 1
                    if (cn > 1.) Cst = (xx*Cst + Csum)/cn
                else                      ! w = negative
                    if (k == ke+1) cycle  ! inflow
                    cn = - dt*w_half(k)/dz(k)
                    kk = k
                    ! extension for Courant numbers > 1
                    if (cn > 1.) then
                        Csum = 0.; dzsum = 0.
                        dtw  = -dt*w_half(k)
                        do while (dzsum+dz(kk) < dtw)
                            if (kk == ks) exit
                            dzsum = dzsum + dz(kk)
                            Csum  =  Csum +  C(kk)
                            kk    =    kk + 1
                        enddo
                        xx = (dtw-dzsum)/dz(kk)
                    else
                        xx = cn
                    endif
                    Cm = C_right(kk) - C_left(kk)
                    C6 = 6.0*(C(kk) - 0.5*(C_right(kk) + C_left(kk)))
                    if (kk == ke) C6 = 0.
                    Cst = C_left(kk) + 0.5*xx*(Cm + (1.0 - tt*xx)*C6)
                    ! extension for Courant numbers > 1
                    if (cn > 1.) Cst = (xx*Cst + Csum)/cn
                endif
                Cst = Cst / dz(kk)  ! for conservation
                flux(k) = w_half(k)*Cst
                ! if (xx > 1.) cflerr = cflerr+1
                ! cflmaxx = max(cflmaxx,xx)
                ! cflmaxc = max(cflmaxc,cn)
                ! }}}
            enddo ! }}}

        case default
            call error_mesg("Not setup diff_method option. &
                             please check input.nml")
    end select

    ! vertical advective tendency
    select case (eqn_form)
        case ("FLUX_FORM")
            do k = ks, ke
                dC_dt     = - ( flux(k+1) - flux(k) ) / dz(k)
                C(k) = C(k) + dC_dt * dt
            end do
        case ("ADVECTIVE_FORM")
            do k = ks, ke
                dC_dt     = - ( flux(k+1) - flux(k) - &
                                C(k)*(w_half(k+1)-w_half(k)) ) / dz(k)
                C(k) = C(k) + dC_dt * dt
            end do
        case default
            call error_mesg("No setup equation form.")
    end select

    end subroutine compute_advection ! }}}


    subroutine slope_z(C, dz, slope, limit, linear) ! {{{
    real, dimension(:), intent(in)  :: C, dz
    real, dimension(:), intent(out) :: slope
    logical,  optional, intent(in)  :: limit, linear

    integer :: k, n
    real    :: grad(2:size(C))
    real    :: Cmin, Cmax
    logical :: limiters = .true.
    logical :: dolinear = .true.

    if (present( limit)) limiters = limit
    if (present(linear)) dolinear = linear

    n = size(C)

    ! compute slope (weighted for unequal levels)
    do k = 2, n
        grad(k) = (C(k)-C(k-1))/(dz(k)+dz(k-1))
    enddo
    if (dolinear) then
        do k = 2, n-1
            slope(k) = (grad(k+1)+grad(k))*dz(k)
        enddo
    else
        do k = 2, n-1
            slope(k) = ( grad(k+1)*(2.*dz(k-1)+dz(k)) + &   ! Equation 1.7
                         grad(k  )*(2.*dz(k+1)+dz(k)) ) * dz(k) &
                     / (   dz(k-1) + dz(k) + dz(k+1)  )
        enddo
    endif
    slope(1) = 2.*grad(2)*dz(1)
    slope(n) = 2.*grad(n)*dz(n)

    ! apply limiters to slope
    if (limiters) then
        do k = 1, n
            if (k >= 2 .and. k <= n-1) then
                Cmin = min(C(k-1), C(k), C(k+1))
                Cmax = max(C(k-1), C(k), C(k+1))
                slope(k) = sign(1.,slope(k)) *  &
                            min( abs(slope(k)), &
                                2.*(C(k)-Cmin), &
                                2.*(Cmax-C(k))  )   ! Equation 1.8
            else
                slope(k) = 0.  ! always slope=0
            endif
        enddo
    endif

    end subroutine slope_z  ! }}}


    subroutine compute_weights(dz, zwt) ! {{{
    real, intent(in),  dimension(:)            :: dz
    real, intent(out), dimension(0:3,size(dz)) :: zwt
    real    :: denom1, denom2, denom3, denom4, num3, num4, x, y
    integer :: k

    do k = 3, size(dz)-1
        denom1 = 1.0/(  dz(k-1) +   dz(k))
        denom2 = 1.0/(  dz(k-2) +   dz(k-1) + dz(k) + dz(k+1))
        denom3 = 1.0/(2*dz(k-1) +   dz(k))  
        denom4 = 1.0/(  dz(k-1) + 2*dz(k))  
        num3   = dz(k-2) + dz(k-1)          
        num4   = dz(k)   + dz(k+1)        
        x      = num3*denom3 - num4*denom4        
        y      = 2.0*dz(k-1)*dz(k)  ! everything up to this point is just
                                    ! needed to compute x1,x1,x3                      
        zwt(0,k) = dz(k-1)*denom1               ! = 1/2 in equally spaced case
        zwt(1,k) = zwt(0,k) + x*y*denom1*denom2 ! = 1/2 in equally spaced case
        zwt(2,k) = dz(k-1)*num3*denom3*denom2   ! = 1/6 ''
        zwt(3,k) = dz(k)*num4*denom4*denom2     ! = 1/6 ''
    enddo

    end subroutine compute_weights  ! }}}

end module dynamics_driver_mod
