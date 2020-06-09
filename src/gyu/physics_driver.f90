module physics_driver_mod
use            global_mod
use          namelist_mod
use          constant_mod
use     error_handler_mod, only: error_mesg
use        substeping_mod, only: time_substeping
use         collision_mod, only: coad1d
contains

    subroutine physics_driver(tidx)
        implicit none
        integer, intent(in) :: tidx

        integer :: k, n
        integer :: num_substep
        real    :: delta_time
        real, dimension(nz)   :: dqv
        real, dimension(nz)   :: dTemp
        real, dimension(nz)   :: qc_
        real, dimension(nbin) :: dm
        real, dimension(nbin) :: Nr_1st_layer

        delta_time = dt
        Nr_1st_layer = Nr(:,1)
        dm = mass_boundary(2:nbin+1,1) - mass_boundary(1:nbin,1)

        vertical_loop: do k = 1, nz 

            ! Calculate mass tendency
            call conc_growth( T(k), qv(k), Prs(k),    & 
                              dm_dt(:,k), dmb_dt(:,k) )

            num_substep = 1
            call time_substeping( dm_dt(:,k), dm, dt, num_substep )

            substeping_loop: do n = 1, num_substep, 1

                ! Compute concentration
                ! TODO: PPM conservation test
                call compute_conc( dm_dt(:,k), dmb_dt(:,k), &
                                   dm, mass(:,k), Nr(:,k) )

            end do substeping_loop

            ! Reinit dt
            dt = delta_time 

            ! Online Coupling with T and qv
            dqv(k)   = - sum( dm_dt(:,k)*dt )
            dTemp(k) = - ( L*dqv(k) ) / ( rho_liquid*Cp )
            qv(k)    = qv(k) + dqv(k)
            T(k)     = T(k)  + dTemp(k)

            ! Change unit for collision
            qc_(k) = sum(Nr(:,k)*mass(:,k)) * rho_liquid ! [kg kg-1] -> [kg m-3]

            ! Compute Stochastic Collision Equation
            call coad1d( dt, nbin, r0, qc_(k), Nr )

            ! TODO: Sedimentation

        end do vertical_loop

        ! Convert T[K] to Theta[K] for dynamics process
        THETA(:) = T(:)*((P0/Prs(:))**(R/Cp))

        ! Reinit distribution for 1st layer
        Nr(:,1) = Nr_1st_layer

    end subroutine physics_driver

    subroutine conc_growth(temp, qv, Pinit, dm_dt, dmb_dt) !{{{
        implicit none
        real,               intent(in)    :: temp, qv, Pinit
        real, dimension(:), intent(out)   :: dm_dt
        real, dimension(:), intent(out)   :: dmb_dt

        real                    :: e, es, RH, S, Fk, Fd
        real, dimension(nbin)   :: Vf
        real, dimension(nbin+1) :: Vfb

        call cal_es_Fk_Fd(temp, Pinit, es, Fk, Fd) 
        e     = Pinit * qv/0.622    ! [hPa] Vapor pressure   
        RH    = (e/es)              ! [%]   Relative humidity
        S     = RH - 1.
        ! S     = 0.01              ! Test

        Vf = 1.; Vfb = 1.
        if (ventilation_effect) then 
            call ventilation(temp, Pinit, radius, Vf)
            call ventilation(temp, Pinit, radius_boundary, Vfb)
        end if

        dm_dt  = 4*PI*radius         *(1./(Fd+Fk))*S*Vf
        dmb_dt = 4*PI*radius_boundary*(1./(Fd+Fk))*S*Vfb

    end subroutine conc_growth  !}}}

    subroutine cal_es_Fk_Fd(temp, Pinit, es, Fk, Fd) !{{{
        implicit none
        real, intent(in)  :: temp, Pinit
        real, intent(out) :: es, Fk, Fd
        
        ! Refer to Rogers & Yau (1996), 16p - (2.17)
        es = 6.112 * exp(( 17.67*(temp-273.15) )/( (temp-273.15)+243.5 ))
        ! To calculate Fd, need to convert the units of 'es'. :: [hPa] > [J m-3]

        ! Refer to Rogers & Yau (1996), 103p - Table 7.1 caption,
        ! "Dv must therefore be multiplied by (1000./P)"        
        Fk = ( (L/(Rv*temp))-1. ) * ( L/(Ka*temp) )    
        Fd = ( Rv*temp ) / ( (Dv*(1000./Pinit)) * (es*100.) )
        
    end subroutine cal_es_Fk_Fd!}}}

    subroutine terminal_velocity(T, P, radius, Vt) !{{{
    ! Input
    ! - T      : Temperature [K]
    ! - P      : Pressure [hPa]
    ! - radius : [m]
    !
    ! Output
    ! - Vt     : Terminal velocity [m s-1]
    ! 
    ! Reference 
    ! - Beard (1976)
        implicit none
        real, intent(in)  :: T, P
        real, intent(in)  :: radius
        real, intent(out) :: Vt

        real :: d0  ! diameter [um]
        real :: T0, P0, l0, mu0
        real :: rho_air, drho
        real :: l, C1, C2, C3, Csc, Cl
        real :: b0, b1, b2, b3, b4, b5, b6
        real :: Da, X, Y, Re, Bo, Np

        d0 = radius * 2 * 1.e6  ! radius [m] -> diameter [um]

        ! constant
        T0  = 293.15     ! [K]
        P0  = 1013.25    ! [hPa]
        l0  = 6.62e-6    ! [cm]
        mu0 = 0.0001818  ! [g cm-1 s-1]

        rho_air    = (P*100.)/(R*T)    ! [kg m-3] air density  ; p = rho R T -> rho = P/RT
        drho = rho_liquid - rho_air    ! drop - air

        ! dynamic viscosity ( Approximate formula, See Yau (1996) - 102p )
        ! mu = 1.72e-5 * ( 393/(T+120.) ) * ( (T/273)**(3./2.) )   ! [kg m-1 s-1]

        l   = l0 * (mu/mu0) * (P0/P) * sqrt(T/T0)  ! [cm] mean free path of air molecules
        C1  = drho * g / (18*mu)        ! [m-1 s-1] = [kg m-3] * [m s-2] / [kg m-1 s-1]
        Csc = 1 + 2.51*l/d0             ! [dimensionless]

        ! Calculate terminal velocity in each regime
        if (d0 < 0.5 ) then
            Vt = 0.     ! ignore
        else if (d0 <= 19) then
            ! Regime 1
            Vt = C1 * Csc * (d0*1.e-6)**2   ! [m s-1] = [m-1 s-1] [dimensionless] [um^2]
        else if (d0 <= 1.07e3) then
            ! Regime 2
            b0 = -0.318657e+1; b1 =  0.992696;    b2 = -0.153193e-2
            b3 = -0.987059e-3; b4 = -0.578878e-3; b5 =  0.855176e-4
            b6 = -0.327815e-5
            C2  = 4 * rho_air * drho* g / ( 3 * mu**2 )  ! [] = [kg2 m-6] [m s-2] / [um2]
            Da  = C2 * (d0*1e-6)**3  ! Davies number [kg2 s-2] = [kg2 m-3 s-2] [m-3]
            X   = log( Da )
            Y   = b0 + b1*X + b2*X**2 + b3*X**3 + b4*X**4 + b5*X**5 + b6*X**6
            Re  = Csc * exp(Y)              ! Reynolds number
            Vt  = mu * Re / (rho_air * d0*1.e-6)
        else !if (d0 <= 7.e3) then
            ! Regime 3
            if (d0 > 7.e3) d0 = 7.e3    ! Assumption: large drop -> d0 = 7.e3 [um]
            b0 = -0.500015e+1; b1 =  0.523778e+1; b2 = -0.204914e+1
            b3 =  0.475294;    b4 = -0.542819e-1; b5 =  0.238449e-2

            ! surface tension [N m-1]
            ! sigma = 7.5 * 1e-2  ; Yau (1996) - 85p
            ! See Yau (1996) - 6.9 problem
            ! Cl = -1.55 * 1e-4   ! [N m-1 K-1]
            ! C2 = 0.118          ! [N m-1]
            ! sigma = Cl*T + C2   ! Resonable at -20 ~ 20 [K] temperature

            C3  = 4 * drho * g / (3*sigma)
            Bo  = C3 * (d0*1.e-6)**2.       ! modified Bond number
            Np  = sigma**3. * rho_air**2. / (mu**4. * drho * g)
            X   = log( Bo * Np**(1./6.) )
            Y   = b0 + b1*X + b2*X**2 + b3*X**3 + b4*X**4 + b5*X**5

            Re  = Np**(1./6.) * exp(Y)
            Vt  = mu * Re / (rho_air*(d0*1.e-6))
        end if  

    end subroutine terminal_velocity!}}}

    subroutine ventilation(T, P, r, Vf) !{{{
    ! Reference
    ! https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
        implicit none
        real,               intent(in ) :: T      ! Temperature [K]
        real,               intent(in ) :: P      ! Pressure    [hPa]
        real, dimension(:), intent(in ) :: r      ! radius      [m]
        real, dimension(:), intent(out) :: Vf     ! ventilation effect

        integer :: i
        real    :: Cd, mu
        real, dimension(size(r)) :: Vt, Re

        ! Assumed constant drag coefficient at large Re.
        ! Cd = 0.45              ! Yau (1996) 125p

        ! Note! We assumed that all drop shape is sphere.
        ! Vt = sqrt( (8./3.)*(r*gravity*rho_liquid)  &
        !                   /(rho_air*Cd) )   ! Yau (1996) equation 8.4
        do i = 1, size(r)
            call terminal_velocity(T, P, r(i), Vt(i))   ! Beard (1976)
        end do

        ! dynamic viscosity of air (See Yau (1996) 102-103p)
        ! mu = 1.72e-5 * ( 393./(T+120.) ) * ( T/273. )**(3./2.)    ! approximate formula

        ! Reynolds number
        Re = 2*rho_air*r*Vt/mu         ! Yau (1996) 116p

        if ( any(0 <= Re .and. Re < 2.5) ) then
            Vf = 1.0 + 0.09*Re
        else ! if (Re > 2.5) then
            Vf = 0.78 + 0.28*sqrt(Re)
        end if

        ! Set maximum value
        where (Vf .gt. 5.0)
            Vf = 5.0
        end where

    end subroutine ventilation !}}}


    subroutine compute_conc( dm_dt, dmb_dt, dm, mass, Nr )   !{{{
        implicit none
        real, dimension(:), intent(in)    :: dm_dt
        real, dimension(:), intent(in)    :: dmb_dt
        real, dimension(:), intent(in)    :: dm
        real, dimension(:), intent(in)    :: mass
        real, dimension(:), intent(inout) :: Nr


        ! TODO: Test physics scheme
        select case (phy_adv_scheme)
            case (1)    ! 1: reassign
                call redistribution( mass, dm_dt, Nr )

            case (2:4)  ! 2: finite_volume, 3: PPM, 4: PPM (but Lin limiter)
                call compute_advection( dt, dmb_dt, dm, Nr )

            case default
                call error_mesg("Not setup phy_adv_scheme option. &
                                 please check input.nml")
        end select

    end subroutine compute_conc !}}}

    subroutine redistribution( mass, dm_dt, Nr ) !{{{
        implicit none
        real, dimension(:) ,intent(in)    :: mass
        real, dimension(:) ,intent(in)    :: dm_dt
        real, dimension(:) ,intent(inout) :: Nr

        integer                   :: i, j, nbin
        real, dimension(size(Nr)) :: next_N
        real, dimension(size(Nr)) :: x, y
        real, dimension(size(mass)) :: next_m

        next_m = mass + dm_dt*dt

        nbin = size(Nr)
        y    = 0.
        ! next_N = 0.

        do i = 1, nbin-1
            ! print*, "i=",i,"mass(i)=",mass(i),"next_m(i)",next_m(i),"mass(i+1)=",mass(i+1)
            if ( next_m(i) <= 0. ) then
                next_N(i) = 0.
            else
                do j = 1, nbin-1
                    if ( mass(j) < next_m(i) ) then
                        if ( mass(j+1) > next_m(i) ) then
                            x(i)      = ( Nr(i)*(next_m(i)-mass(j+1)) ) &
                                      / (          mass(j)-mass(j+1) )
                            y(i+1)    = Nr(i) - x(i)
                            next_N(i) = x(i) + y(i)
                            exit
                        end if
                    end if
                end do
            end if
            ! print*,"i=",i,"Nr(i)=",Nr(i),"x=",x(i),"y=",y(i),"next_N(i)=",next_N(i)
            ! print*, mass(i),mass(i+1),next_N(i)
        end do
        next_N(nbin) = Nr(nbin)+y(nbin)
        
        ! print*, "TODO! redistribution"
        ! print*, "     Porting from ncl code... (+ doc/redistribution/reassign.ncl)"
        ! stop

        Nr = next_N

    end subroutine redistribution   !}}}

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
    real,               intent(in   ) :: dt
    real, dimension(:), intent(in   ) :: w_half
    real, dimension(:), intent(in   ) :: dz
    real, dimension(:), intent(inout) :: C

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
    ! if (phy_adv_scheme == 1) do_outflow_bnd = .false.

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

    select case (phy_adv_scheme)
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

                ! TODO: Test courant error while substeping ...
                ! if (cn > 1.) call error_mesg("Courant number > 1 &
                !                               in physics process")
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
                flux(k) = w_half(k)*Cst
                ! if (xx > 1.) cflerr = cflerr+1
                ! cflmaxx = max(cflmaxx,xx)
                ! cflmaxc = max(cflmaxc,cn)
                ! }}}
            enddo ! }}}

        case default
            call error_mesg("Not setup physics_driver_mod option. &
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
            call error_mesg("No setup physics equation form.")
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

end module physics_driver_mod
