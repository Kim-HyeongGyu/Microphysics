module microphysics_mod
use            global_mod !only: error_mesg
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
        mass(:,1)          = m
        mass_boundary(:,1) = mb

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

    subroutine compute_mass()
        ! 1) interpolation
        ! 2) advection
        select case (mass_scheme)
            case ("interpolation")
                stop
            case ("FVM")
                ! call mass_advection(dmb_dt, N, dt, nbin, &
                !                   dm, "finite_volume", next_N)
                stop
            case ("PPM")
                ! call mass_advection(dmb_dt, N, dt, nbin, &
                !                   dm, "finite_volume", next_N)
                stop
            case default
                
        end select 
    end subroutine compute_mass

end module microphysics_mod
