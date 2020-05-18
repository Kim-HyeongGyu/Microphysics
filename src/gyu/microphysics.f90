module microphysics_mod
use            global_mod !only: error_mesg
use       conc_advect_mod, only: conc_advection
contains
    subroutine conc_dist() !{{{
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

    end subroutine conc_dist!}}}

    subroutine make_bins() !{{{
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

        allocate(radius(nbin), radius_boundary(nbin+1))
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

    end subroutine make_bins!}}}

    subroutine conc_growth(temp, qv, Pinit, dm_dt, dmb_dt) !{{{
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
            call ventilation(temp, Pinit, radius, Vf)
            call ventilation(temp, Pinit, radius_boundary, Vfb)
        end if
        dm_dt  = 4*PI*radius*(1./(Fd+Fk))*S*Vf
        dmb_dt = 4*PI*radius_boundary*(1./(Fd+Fk))*S*Vfb
        !do i = 1, nbin
        !print*,dm_dt(i)
        !end do

    end subroutine conc_growth!}}}

    subroutine cal_es_Fk_Fd(temp, Pinit, es, Fk, Fd) !{{{
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
        real :: R, g, T0, P0, l0, mu0
        real :: rho_liquid, rho_air, drho, mu
        real :: l, C1, C2, C3, Csc, Cl
        real :: b0, b1, b2, b3, b4, b5, b6
        real :: Da, X, Y, Re, sigma, Bo, Np

        d0 = radius * 2 * 1.e6  ! radius [m] -> diameter [um]

        ! constant
        R   = 287.       ! [J kg-1 K-1]
        g   = 9.8        ! [m s-2]
        T0  = 293.15     ! [K]
        P0  = 1013.25    ! [hPa]
        l0  = 6.62e-6    ! [cm]
        mu0 = 0.0001818  ! [g cm-1 s-1]

        rho_liquid = 1000.      ! [kg m-3] water density
        rho_air    = (P*100.)/(R*T)    ! [kg m-3] air density  ; p = rho R T -> rho = P/RT
        drho = rho_liquid - rho_air    ! drop - air

        ! dynamic viscosity ( Approximate formula, See Yau (1996) - 102p )
        mu = 1.72e-5 * ( 393/(T+120.) ) * ( (T/273)**(3./2.) )   ! [kg m-1 s-1]

        l   = l0 * (mu/mu0) * (P0/P) * sqrt(T/T0)  ! [cm] mean free path of air molecules
        C1  = drho * g / (18*mu)        ! [m-1 s-1] = [kg m-3] * [m s-2] / [kg m-1 s-1]
        Csc = 1 + 2.51*l/d0             ! [dimensionless]

        ! Calculate terminal velocity in each regime
        if (d0 .lt. 0.5 ) then
            Vt = 0.     ! ignore
        else if (d0 .le. 19) then
            ! Regime 1
            Vt = C1 * Csc * (d0*1.e-6)**2   ! [m s-1] = [m-1 s-1] [dimensionless] [um^2]
        else if (d0 .le. 1.07e3) then
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
        else if (d0 .le. 7.e3) then
            ! Regime 3
            b0 = -0.500015e+1; b1 =  0.523778e+1; b2 = -0.204914e+1
            b3 =  0.475294;    b4 = -0.542819e-1; b5 =  0.238449e-2

            ! surface tension [N m-1]
            ! sigma = 7.5 * 1e-2  ; Yau (1996) - 85p
            ! See Yau (1996) - 6.9 problem
            Cl = -1.55 * 1e-4   ! [N m-1 K-1]
            C2 = 0.118          ! [N m-1]
            sigma = Cl*T + C2   ! Resonable at -20 ~ 20 [K] temperature

            C3  = 4 * drho * g / (3*sigma)
            Bo  = C3 * (d0*1.e-6)**2.       ! modified Bond number
            Np  = sigma**3. * rho_air**2. / (mu**4. * drho * g)
            X   = log( Bo * Np**(1./6.) )
            Y   = b0 + b1*X + b2*X**2 + b3*X**3 + b4*X**4 + b5*X**5

            Re  = Np**(1./6.) * exp(Y)
            Vt  = mu * Re / (rho_air*(d0*1.e-6))
        else
            ! terminal velocity [m s-1] at d0 = 7.e3 [um] (using regime 3)
            d0 = 7.e3
            b0 = -0.500015e+1; b1 =  0.523778e+1; b2 = -0.204914e+1
            b3 =  0.475294;    b4 = -0.542819e-1; b5 =  0.238449e-2

            Cl = -1.55 * 1e-4   ! [N m-1 K-1]
            C2 = 0.118          ! [N m-1]
            sigma = Cl*T + C2   ! Resonable at -20 ~ 20 [K] temperature

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
        real, dimension(size(r)) :: Vt, Re
        real :: rho_liquid, rho_air, gravity
        real :: Cd, mu

        rho_liquid = 1000.     ! [kg m-3] water density
        rho_air    = 1.225     ! [kg m-3] air density
        gravity    = 9.8

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
        mu = 1.717e-5          ! [kg m-1 s-1] (at 273 [K])

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

    subroutine redistribution(Nr, next_N, mass, next_m) !{{{
        implicit none
        real, dimension(nbin) :: Nr, mass, next_m, next_N
        real, dimension(nbin) :: x, y
        integer :: i

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

    end subroutine redistribution!}}}

end module microphysics_mod
