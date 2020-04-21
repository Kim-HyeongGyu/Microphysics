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
        ! TODO: add level index
        mass(:,1)          = m
        mass_boundary(:,1) = mb

    end subroutine make_bins

    subroutine conc_growth(temp, qv, Pinit, dmdt)
        implicit none
        real, dimension(nz), intent(in)  :: temp, qv, Pinit
        real, dimension(nz), intent(out) :: dmdt

        real :: Rv, Dv, Ka, L
        real, dimension(nz) :: e, es, RH, S, Fd, Fk

        ! S     = RH - 1.
        S     = 0.01                ! For test
        ! temp  = 293.15              ! For test             [K]
        Rv    = 462                 ! vapor gas constant   [J kg-1 K-1]
        L     = 2.5e6               ! heat of vaporization [W m-2]
        ! Refer to Yau (1996), 103p - Table 7.1
        Dv    = 2.21e-5             ! Diffusion            [m2 s-1]
        Ka    = 2.40e-2             ! Thermal conductivity [J m-1 s-1 K-1]

        ! Refer to 대기열역학
        e     = Pinit * qv/0.622    ! vapor pressure       [hPa]
        es    = 6.112 * exp(( 17.67*temp )/( temp+243.5 )) ! saturated e
        Rh    = (e/es)*100          ! Relative humidity    [%]

        Fd    = ( Rv*temp ) / (Dv*es)  
        Fk    = ( L/(Rv*temp) - 1. ) * ( L/(Ka*temp) )
        dmdt  = 4*PI*radius*(1./(Fd+Fk))*S
        
    end subroutine conc_growth

end module microphysics_mod
