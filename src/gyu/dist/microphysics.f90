module microphysics
contains
    subroutine conc_dist(radius, Nc, qc, dist_type, Nr)
        real, dimension(:), intent(in) :: radius   ! [um]
        real,               intent(in) :: Nc, qc
        character(len=*),   intent(in) :: dist_type
        real, dimension(size(radius)), intent(out) :: Nr
        real :: u, std, lamda, N0, r0
        real :: umul
        real, parameter :: rho = 1000.     ! [kg m-3] water density
        real, parameter :: PI = 3.141592

        select case (trim(dist_type))
            ! 1) Log normal dist_type
            case ("log_normal")
                N0  = Nc
                r0  = ( (3./(4.*PI))*(qc/Nc)*(1./rho) )**(1./3.)
                u   = log(r0)
                std = 2.5
                ! TODO! std = 0 when use below formula
                ! std = sqrt((1./3.)* &
                !           ((2./3.)*log((3./4.)*qc/(PI*rho*N0)) -2*u))
                Nr  = N0 / (sqrt(2*PI)*std*radius) &
                    * exp(-( (log(radius)-u)**2 / (2*std**2) ))
            ! 2) Gamma dist_type
            case ("gamma")
                ! u     = 1.e9/Nc   ! Note small Nc makes infinity error
                u     = 10.         ! usually, 2~15
                umul  = (u+1)*(u+2)*(u+3)
                N0    = ( Nc/gamma(u+1) ) &
                      * ( (rho/qc)*(4./3.)*PI*Nc*umul )**((u+1)/3.)
                lamda = ( (rho/qc)*(4./3.)*PI*Nc*umul )**(   1./3.)
                Nr    = N0*radius**u*exp(-lamda*radius)
            case default
                ! call error_mesg("Not setup dist_type option. &
                !                  please check input.nml")
                stop 2
        end select
    end subroutine conc_dist

    subroutine make_bins(rmin, rratio, nbin, radius)
        integer,               intent(in)  :: nbin
        real,                  intent(in)  :: rmin, rratio
        real, dimension(nbin), intent(out) :: radius ! [m]
        real, dimension(nbin)   :: r    ! radius [m]
        real, dimension(nbin)   :: m    ! mass   [kg]
        real, dimension(nbin+1) :: rb   ! radius at boundary
        real, dimension(nbin+1) :: mb   ! mass at boundary
        real, parameter :: PI = 3.141592
        real, parameter :: rho = 1000.  ! [kg m-3] water density
        
        rb = (/ (rmin*(rratio**i), i=0,nbin) /)
        mb = (4./3.)*pi*rho*rb**3

        ! Interpolate using mass
        do i = 1, nbin
            m(i) = ( mb(i)+mb(i+1) ) / 2.
            r(i) = ( (3./4)/(PI*rho)*m(i) )**(1./3.)
        enddo
        radius = r
    end subroutine make_bins

end module microphysics
