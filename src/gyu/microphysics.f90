module microphysics_mod
use          constant_mod
use          namelist_mod
use     error_handler_mod, only: error_mesg
contains
    subroutine initialize_microphysics( radius, radius_boundary, &
                                        mass, mass_boundary,     & 
                                        Nr, dm_dt, dmb_dt, r0 )
        real, dimension(  :), allocatable, intent(out) :: radius         
        real, dimension(  :), allocatable, intent(out) :: radius_boundary
        real, dimension(:,:), allocatable, intent(out) :: mass           
        real, dimension(:,:), allocatable, intent(out) :: mass_boundary  
        real, dimension(:,:), allocatable, intent(out) :: Nr             
        real, dimension(:,:), allocatable, intent(out) :: dm_dt
        real, dimension(:,:), allocatable, intent(out) :: dmb_dt
        real,                              intent(out) :: r0

        call make_bins(radius, radius_boundary, mass, mass_boundary)
        call conc_dist(Nc, qc, radius, radius_boundary, Nr, r0)

        allocate( dm_dt(nbin  ,nz))
        allocate(dmb_dt(nbin+1,nz))

    end subroutine initialize_microphysics

    subroutine make_bins(radius, radius_boundary, mass, mass_boundary) !{{{
        real, dimension(  :), allocatable, intent(out) :: radius    ! [m]
        real, dimension(  :), allocatable, intent(out) :: radius_boundary
        real, dimension(:,:), allocatable, intent(out) :: mass      ! [kg]
        real, dimension(:,:), allocatable, intent(out) :: mass_boundary

        real, dimension(nbin  ) :: m    ! [kg] mass   
        real, dimension(nbin+1) :: mb   ! mass at boundary

        if ( drop_var == 1 ) then
            print*, "::: droplet variable used: rmin & rmax :::"
            rratio = (rmax/rmin)**(1./nbin)
        elseif ( drop_var == 2 ) then
            print*, "::: droplet variable used: rmin & rratio :::"
        endif

        allocate(         radius(nbin  )   )
        allocate(radius_boundary(nbin+1)   )
        allocate(           mass(  nbin,nz))
        allocate(  mass_boundary(nbin+1,nz))

        radius_boundary = (/ (rmin*(rratio**i), i=0,nbin) /)
        mb = (4./3.)*pi*rho_liquid*radius_boundary**3

        ! Interpolate using mass
        do i = 1, nbin
            m(i) = ( mb(i)+mb(i+1) ) / 2.
            radius(i) = ( (3./4)/(PI*rho_liquid)*m(i) )**(1./3.)
        enddo

        do k = 1, nz 
            mass(:,k)          = m
            mass_boundary(:,k) = mb
        end do

    end subroutine make_bins    !}}}

    subroutine conc_dist(Nc, qc, radius, radius_boundary, Nr2d, r0) !{{{
        implicit none
        real,                              intent( in) :: Nc, qc
        real, dimension(:  ),              intent( in) :: radius
        real, dimension(:  ),              intent( in) :: radius_boundary
        real, dimension(:,:), allocatable, intent(out) :: Nr2d
        real,                              intent(out) :: r0

        integer :: i, k
        real :: u, std, lambda, N0
        real :: umul, r0_max
        real, dimension(nbin) :: Nr ! Initial number of droplet
        real, dimension(nbin) :: m  ! mass
        real, dimension(nbin) :: dr
        real, dimension(nbin) :: dN_dlnr

        allocate(Nr2d(nbin,nz))

        dr = radius_boundary(2:) - radius_boundary(:nbin-1)

        select case (dist_type)
            case (1)    ! 1: Log normal distribution
                N0  = Nc
                r0  = 1.e-5
                r0_max = ( qc / (Nc*rho_liquid*(4./3.)*pi) )**(1./3.)
                if (r0 > r0_max) then
                    ! print*, "Note! r0 > r0_max"
                    ! print*, "r0: ", r0, ",  ", "r0_max: ", r0_max
                    r0 = r0_max*0.9
                    ! print*, "r0 will be changed to => ", r0
                end if    

                u   = log(r0)
                std = sqrt((2./9.)*log(qc/(rho_liquid*(4./3.)*pi*(r0**3.)*Nc)))
                Nr  = N0/(sqrt(2*pi)*radius*std)      &
                    * exp(- ( (log(radius)-u)**2 )    &
                          / (    2*std**2 ) ) * dr
                m   = (4./3.)*pi*rho_liquid*(radius**3)*Nr

            case (2)    ! 2: Gamma distribution
                u      = min(1.e9/Nc+2., 15.)  ! usually, 2~15
                umul   = (u+1)*(u+2)*(u+3)
                lambda = ( (Nc/qc)*(4./3.)*PI*umul*rho_liquid )**(1./3.)
                Nr     = Nc / gamma(u+1)                &
                            * lambda*(lambda*radius)**u &
                            * exp(-lambda*radius)*dr
                m      = Nr*rho_liquid*pi*(4./3.)*radius**3.
                r0     = u / lambda

            case default
                call error_mesg("Not setup dist_type option. &
                                 please check input.nml")
        end select

        do i = 1, nbin
            dN_dlnr(i) = Nr(i)/log( radius_boundary(i+1) &
                                  / radius_boundary(i)   )
        end do

        ! do k = 1, nz
        !     Nr2d(:,k) = Nr
        ! end do
        Nr2d = 0.
        Nr2d(:,1) = Nr  ! Give distribution at 1st layer

    end subroutine conc_dist!}}}

end module microphysics_mod
