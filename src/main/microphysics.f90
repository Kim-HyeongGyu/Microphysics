module microphysics_mod
use            global_mod !only: error_mesg
use         advection_mod, only: slope_z, compute_weights
contains
    subroutine conc_dist()
        real :: u, std, lambda, N0, r0
        real :: umul, r0_max
        real, dimension(nbin) :: m  ! mass
        real, dimension(nbin) :: dr

        allocate(Nr(nbin), dN_dlnr(nbin))

        dr = radius_boundary(2:) - radius_boundary(:nbin-1)

        select case (trim(dist_type))
            ! 1) Log normal dist_type {{{
            case ("log_normal")
                N0  = Nc
                r0  = 1.e-5
                r0_max = ( qc / (Nc*rho*(4./3.)*pi) )**(1./3.)
                if (r0 > r0_max) then
                    r0 = r0_max/10.
                    print*,"Note! r0 > r0_max"
                    print*,"r0 will be changed to => ", r0
                end if    

                u   = log(r0)
                std = sqrt((2./9.)*log(qc/(rho*(4./3.)*pi*(r0**3.)*Nc)))
                Nr  = N0/(sqrt(2*pi)*radius*std)      &
                    * exp(- ( (log(radius)-u)**2 )    &
                          / (    2*std**2 ) ) * dr
                m   = (4./3.)*pi*rho*(radius**3)*Nr     !}}}  

            ! 2) Gamma dist_type {{{
            case ("gamma_dist")
                u      = min(1.e9/Nc+2., 15.)  ! usually, 2~15
                umul   = (u+1)*(u+2)*(u+3)
                lambda = ( (Nc/qc)*(4./3.)*PI*umul*rho )**(1./3.)
                Nr     = Nc / gamma(u+1)                &
                            * lambda*(lambda*radius)**u &
                            * exp(-lambda*radius)*dr
                m      = Nr*rho*pi*(4./3.)*radius**3.    !}}}

            case default
                call error_mesg("Not setup dist_type option. &
                                 please check input.nml")
        end select

        do i = 1, nbin
            dN_dlnr(i) = Nr(i)/log( radius_boundary(i+1) &
                                  / radius_boundary(i)   )
        end do

        ! print*, "radius = ", radius
        ! print*, "Nc     = ", sum(Nr)
        ! print*, "qc     = ", sum(m)

    end subroutine conc_dist

    subroutine make_bins()
        real, dimension(nbin)   :: r    ! radius [m]
        real, dimension(nbin)   :: m    ! mass   [kg]
        real, dimension(nbin+1) :: rb   ! radius at boundary
        real, dimension(nbin+1) :: mb   ! mass at boundary

        if ( drop_var == 1 ) then
            print*, "::: droplet variable used: rmin & rmax :::"
            rratio = (rmax/rmin)**(1./nbin)
        elseif ( drop_var == 2 ) then
            print*, "::: droplet variable used: rmin & rratio :::"
        endif

        allocate(radius(nbin), radius_boundary(nbin))
        rb = (/ (rmin*(rratio**i), i=0,nbin) /)
        mb = (4./3.)*pi*rho*rb**3

        ! Interpolate using mass
        do i = 1, nbin
            m(i) = ( mb(i)+mb(i+1) ) / 2.
            r(i) = ( (3./4)/(PI*rho)*m(i) )**(1./3.)
        enddo

        radius             = r
        radius_boundary    = rb
        do k = 1, nz 
            mass(:,k,1)          = m
            mass_boundary(:,k,1) = mb
        end do

    end subroutine make_bins

    subroutine conc_growth(temp, qv, Pinit, dm_dt, dmb_dt)
        implicit none
        real,                    intent(in)  :: temp, qv, Pinit
        real, dimension(nbin),   intent(out) :: dm_dt
        real, dimension(nbin+1), intent(out) :: dmb_dt

        real                    :: e, es, RH, S, Fk, Fd
        real, dimension(nbin)   :: Vf
        real, dimension(nbin+1) :: Vfb

        call cal_es_Fk_Fd(temp,Pinit,es,Fk,Fd) 
        e     = Pinit * qv/0.622    ! vapor pressure       [hPa]
        RH    = (e/es)*100          ! Relative humidity    [%]
        
        ! S     = RH - 1.
        S     = 0.01                ! For test

        Vf = 1.; Vfb = 1.
        if (ventilation_effect) then 
            call ventilation(temp, radius, Vf)
            call ventilation(temp, radius_boundary, Vfb)
        end if
        dm_dt  = 4*PI*radius*(1./(Fd+Fk))*S*Vf
        dmb_dt = 4*PI*radius_boundary*(1./(Fd+Fk))*S*Vfb
        !do i = 1, nbin
        !print*,dm_dt(i)
        !end do

    end subroutine conc_growth

    subroutine cal_es_Fk_Fd(temp, Pinit, es, Fk, Fd)
        implicit none
        real, intent(in)  :: temp, Pinit
        real, intent(out) :: es, Fk, Fd
        
        real :: L, Rv, Ka, Dv
        
        L  = 2500297.8  ! heat of vaporization [J kg-1]
        Rv = 467        ! vapor gas constant   [J kg-1 K-1]

        ! Refer to Rogers & Yau (1996), 103p - Table 7.1 (T=273K)
        Ka = 2.40e-2    ! coefficient of thermal conductivity of air    [J m-1 s-1 K-1]
        Dv = 2.21e-5    ! molecular diffusion coefficient               [m2 s-1]

        es = 6.112 * exp(( 17.67*(temp-273.15) )/( (temp-273.15)+243.5 ))
        ! To calculate Fd, need to convert the units of 'es'. :: [hPa] > [J m-3]

        Fk = ((L/(Rv*temp))-1.)*((L*rho)/(Ka*temp))
        Fd = (rho*Rv*temp) / ((Dv*(1000./Pinit))*(es*100.))
        
    end subroutine cal_es_Fk_Fd

    subroutine ventilation(T, r, Vf)
    ! Reference
    ! https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
        implicit none
        real,               intent(in   ) :: T      ! Temperature [K]
        real, dimension(:), intent(in   ) :: r      ! radius      [m]
        real, dimension(:), intent(inout) :: Vf     ! ventilation effect
        real, dimension(size(r)) :: Vt, Re
        real :: rho_liquid, rho_air, gravity
        real :: Cd, mu

        rho_liquid = 1000.     ! [kg m-3] water density
        rho_air    = 1.225     ! [kg m-3] air density
        gravity    = 9.8

        ! Assumed constant drag coefficient at large Re.
        Cd = 0.45              ! Yau (1996) 125p

        ! Note! We assumed that all drop shape is sphere.
        Vt = sqrt( (8./3.)*(r*gravity*rho_liquid)  &
                          /(rho_air*Cd) )   ! Yau (1996) equation 8.4

        ! dynamic viscosity of air (See Yau (1996) 102-103p)
        mu = 1.72e-5 * ( 393./(T+120.) ) * ( T/273. )**(3./2.)    ! approximate formula
        ! mu = 1.717e-5          ! [kg m-1 s-1] (at 273 [K])

        ! Reynolds number
        Re = 2*rho_air*r*Vt/mu         ! Yau (1996) 116p

        if ( any(0 <= Re .and. Re < 2.5) ) then
            Vf = 1.0 + 0.09*Re
        else ! if (Re > 2.5) then
            Vf = 0.78 + 0.28*sqrt(Re)
        end if
        
    end subroutine ventilation

    subroutine compute_conc(dmb_dt, Nr, next_Nr, mass, next_m) !{{{
    implicit none
    real, dimension(nbin+1), intent(in)  :: dmb_dt
    real, dimension(nbin),   intent(in)  :: Nr, mass, next_m
    real, dimension(nbin),   intent(out) :: next_Nr
    real, dimension(nbin) :: dm
        ! 1) reassign (redistribution)
        ! 2) advection
        dm = mass_boundary(2:nbin+1,1,1) - mass_boundary(1:nbin,1,1)
        select case (mass_scheme)
            case ("reassign")
                call redistribution(Nr, next_Nr, mass, next_m)
            case ("finite_volume","PPM")
                ! Conserve number
                call conc_advection(dmb_dt, Nr, dt, nbin, &
                                        dm, mass_scheme, next_Nr)
            case default
                call error_mesg("Not setup mass_scheme option. &
                                 please check input.nml")
        end select

    end subroutine compute_conc !}}}

    subroutine conc_advection(w_half, C, dt, nz, &   ! {{{
                                 dz, scheme, next_C)

    implicit none
    integer,              intent(in) :: dt, nz
    character(len=*),     intent(in) :: scheme
    real, dimension(nz),  intent(in) :: C, dz
    real, dimension(nz), intent(out) :: next_C
    integer :: k, kk
    integer :: ks, ke, kstart, kend
    real    :: wgt   ! weight dz
    real    :: Cwgt  ! weight variable
    real    :: dC_dt, zbottom = 0.
    real    :: tt, cn, Csum, dzsum, dtw
    real    :: xx, a, b, Cm, C6, Cst, Cdt
    real, dimension(0:3,nz) :: zwt
    real, dimension(nz)    :: slp, C_left, C_right
    real, dimension(nz)   :: slope
    real, dimension(nz+1) :: w_half, flux
    character(len=20)     :: eqn_form = "FLUX_FORM"
    logical :: linear, test_1

    ! vertical indexing
    ks     =    1; ke   = nz
    ! kstart = ks+1; kend = ke
    kstart =   ks; kend = ke+1  ! do outflow boundary

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
                    Cst = Cst / dz(k-1)
                else
                    if (k == ke+1) cycle        ! inflow
                    cn  = -dt*w_half(k)/dz(k)
                    Cst = C(k) - 0.5*slope(k)*(1.-cn)
                    Cst = Cst / dz(k)
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
                Cst = Cst / dz(kk)  ! for conservation
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
                dC_dt     = - ( flux(k+1) - flux(k) )
                next_C(k) = C(k) + dC_dt * dt
! if (next_C(k) >= 1.e100) then
!     print*, k, next_C(k), C(k), dC_dt, dt
! end if
            end do
        case ("ADVECTIVE_FORM")
            do k = ks, ke
                ! dC_dt     = - ( flux(k+1)/dz(k+1) - flux(k)/dz(k) ) &
                !             - ( C(k)*(w_half(k+1)-w_half(k)) ) / dz(k)
                dC_dt     = - ( flux(k+1) - flux(k) ) &
                            - ( C(k)*(w_half(k+1)-w_half(k)) ) / dz(k)
                next_C(k) = C(k) + dC_dt * dt
            end do
        case default
            call error_mesg("No setup equation form.")
    end select  ! }}}

    end subroutine conc_advection ! }}}

    subroutine redistribution(Nr, next_N, mass, next_m)
        implicit none
        real, dimension(nbin) :: Nr, mass, next_m, next_N
        real, dimension(nbin) :: x, y

        y = 0
        do i = 1, nbin-1
        ! print*, "i=",i,"mass(i)=",mass(i),"next_m(i)",next_m(i),"mass(i+1)=",mass(i+1)
            if (mass(i) < next_m(i)) then
                if(mass(i+1) > next_m(i)) then
                    x(i)=(Nr(i)*(next_m(i)-mass(i+1)))/(mass(i)-mass(i+1))
                    y(i)=Nr(i)-x(i)
                    if (i==1) then
                        next_N(i) = x(i)
                    else
                        next_N(i)=x(i)+y(i-1)
                    end if
                else
                    Nr(i+1) = Nr(i+1)+Nr(i)
                    next_N(i)=0
                end if
            else
                next_N(i)=Nr(i)
            end if
            ! print*,"i=",i,"Nr(i)=",Nr(i),"Nr(i+1)=",Nr(i+1),"next_N(i)=",next_N(i)
        ! print*, mass(i),mass(i+1),next_N(i)
        end do

        next_N(nbin) = Nr(nbin)+y(nbin-1)
        !print*,Nr(nbin),y(nbin-1),next_N(nbin)
        ! Nr = next_N
        !print*, Nr
        !print*, sum(Nr)

    end subroutine redistribution

end module microphysics_mod
