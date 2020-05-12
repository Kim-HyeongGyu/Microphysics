module advection_mod
use            global_mod !only: error_mesg
contains
    subroutine compute_advection(w_full, C, dt, nz, &   ! {{{
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
    integer :: kk
    integer :: ks, ke, kstart, kend
    real    :: wgt   ! weight dz
    real    :: Cwgt  ! weight variable
    real    :: dC_dt, zbottom = 0.
    real    :: tt, cn, Csum, dzsum, dtw
    real    :: xx, a, b, Cm, C6, Cst, Cdt
    real, dimension(0:3,nz) :: zwt
    real, dimension(nz)    :: slp, C_left, C_right
    real, dimension(nz)   :: w_full, slope
    real, dimension(nz+1) :: w_half, flux
    character(len=20)     :: eqn_form = "FLUX_FORM"
    logical :: linear, test_1

    ! vertical indexing
    ks     =    1; ke   = nz
    ! kstart = ks+1; kend = ke
    kstart =   ks; kend = ke+1  ! do outflow boundary

    ! Make stagged grid for advection
    do k = ks+1, ke
        wgt = dz(k-1) / ( dz(k-1)+dz(k) )
        w_half(k) = w_full(k-1) + wgt*(w_full(k)-w_full(k-1))
    end do

    ! TODO: How to give boundary condition?
    ! Set Boundary Condition
    ! most likely w = 0 at these points
     w_half(1) = 0.; w_half(nz+1) = 0.     ! Homogeneous Dirichlet BC
    ! w_half(1) = 2.; w_half(nz+1) = 2.
    ! flux(ks) = 0.; flux(ke+1) = 0.        ! Neumann BC
    flux(ks)   = w_half(ks)*C(ks)           ! do outflow boundary
    flux(ke+1) = w_half(ke+1)*C(ke)

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

                if (cn > 1.) call error_mesg("Courant number > 1")
            end do !}}}

        ! 3) Piecewise Parabolic Method, Colella and Woodward (1984) {{{
        case ("PPM")
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

            ! if (diff_scheme == FINITE_VOLUME_PARABOLIC2) then
            ! limiters from Lin (2003), Equation 6 (relaxed constraint)
                ! do k = ks, ke
                ! C_left (k) = C(k) - sign( min(abs(slp(k)),       &
                !                           abs(C_left (k)-C(k))), &
                !                           slp(k) )
                ! C_right(k) = C(k) + sign( min(abs(slp(k)),       &
                !                           abs(C_right(k)-C(k))), &
                !                           slp(k) )
                ! enddo
            ! else
            ! limiters from Colella and Woodward (1984), Equation 1.10
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
            ! endif
            
            ! compute fluxes at interfaces
            tt = 2./3.
            do k = kstart, kend ! ks+1, nz
                if (w_half(k) >= 0.) then ! w = positive {{{
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
                    if (cn > 1.) Cst = (xx*Cst + Csum)/cn   ! }}}
                else                        ! w = negative {{{
                    if (k == ke+1) cycle    ! inflow
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
                endif   ! }}}
                flux(k) = w_half(k)*Cst
                ! if (xx > 1.) cflerr = cflerr+1
                ! cflmaxx = max(cflmaxx,xx)
                ! cflmaxc = max(cflmaxc,cn)
            enddo
            ! }}}

        case default
            call error_mesg("Not setup diff_method option. &
                             please check input.nml")
    end select
    
    ! vertical advective tendency {{{
    select case (eqn_form)
        case ("FLUX_FORM")
            do k = ks, ke
                ! Note: for conserve quantity, dz index is different.
                !      Discuss with minwoo (2020.04.20)
                dC_dt     = - ( flux(k+1)/dz(k+1) - flux(k)/dz(k) )
                next_C(k) = C(k) + dC_dt * dt
            end do
        case ("ADVECTIVE_FORM")
            do k = ks, ke
                dC_dt     = - ( flux(k+1)/dz(k+1) - flux(k)/dz(k) ) &
                            - ( C(k)*(w_half(k+1)-w_half(k)) ) / dz(k)
                next_C(k) = C(k) + dC_dt * dt
            end do
        case default
            call error_mesg("No setup equation form.")
    end select  ! }}}

    end subroutine compute_advection ! }}}


    subroutine slope_z(C, dz, slope, limit, linear) ! {{{
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
                                2.*(Cmax-C(k))  )   ! Equation 1.8
            else
                slope(k) = 0.  ! always slope=0
            endif
        enddo
    endif

    end subroutine slope_z  ! }}}


    subroutine compute_weights(dz, zwt) ! {{{
    real, intent(in),  dimension(:)      :: dz
    real, intent(out), dimension(0:3,nz) :: zwt
    real    :: denom1, denom2, denom3, denom4, num3, num4, x, y
    integer :: k
   
    do k = 3, nz-1
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

end module advection_mod
