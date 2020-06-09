module collision_mod
! -----------------------------------------------------------
!
! stochastic collection equation (
! code :
! https://www2.meteo.uni-bonn.de/forschung/gruppen/tgwww/people/abott/fortran/fortran_english.html
! Reference: Bott, A., 1998.
!
! -----------------------------------------------------------
implicit none

contains

subroutine coad1d( dt, n, rq0_, xmw_, nr )
    implicit none
    real,               intent(in)  :: dt
    integer,            intent(in)  :: n
    real,               intent(in)  :: rq0_, xmw_
    real, dimension(n), intent(out) :: nr
     
    ! integer, parameter :: n = 400   <- nbin
    real                 :: rq0, xmw
    real                 :: dlnr, ax, xn0, xn1 ,x1, &
                                       x0, lmin, t, tlmin
    real, parameter      :: scal = 2.
    real, parameter      :: emin = 1.e-9, pi = 3.141592654
    real, dimension(n,n) :: c, ima, ck, ec
    real, dimension(n)   :: g, r, e
    integer :: i, imax, j
    integer :: isw = 1  ! Hall (1980) kernel 
    !   implicit double precision (a-h,o-z)
    !common /grid/ g(n),r(n),e(n)
    !common /cour/ c(n,n),ima(n,n)
    !common /kern/ ck(n,n),ec(n,n)
    !dimension rri(n),eei(n)
    !data emin,pi,tmax/1.e-9,3.141592654,3600./

! dt  : time step (input in sec)
! g   : spectral mass distribution (mg/cm**3)
! e   : droplet mass grid (mg)
! r   : droplet radius grid (um)
! dlnr: constant grid distance of logarithmic grid 
! rq0 : mode radius of initial distribution (input in um)  <- change [m]
! xmw : total water content (input in g/m**3)              <- change [kg m-3]
! xn0 : mean initial droplet mass (mg)
! xn1 : total initial droplet number concentration (1/cm^3)
! ax  : growth factor for consecutive masses
! scal: scaling factor for calculation of ax, see below
! isw : collision kernel: 0 (long), 1 (hall), 2 (golovin)  <- Hall kernel used
! read input variables
    ! Mean drop radius (See 2287p)                  !      _                            
    rq0  = rq0_*1.e+02                              ! rq0: r  [m] -> [cm]  
    ! Toal cloud water content (See 2287p) 
    xmw  = xmw_*1.e-6                               ! xmw: L [kg m-3] -> [mg cm-3]
    ! Mean drop mass (See 2287p)                    !      _
    xn0  = (4./3.)*pi*1000.*exp(log(rq0)*3.)        ! xn0: x _ [mg]
    xn1  = xmw/xn0                                  ! xn1: L/x [cm-3]              
    dlnr = log(2.e0)/(3.*scal)                      ! dlnr: dy = ln(alpha/3)
    ax   = 2.e0**(1.0/scal)                         ! ax:   alpha = 2^(1/2) 
    
    ! mass and radius grid
    e(1) = emin*0.5*(ax+1.)                         ! e: Minimum mass x1=x/2
    r(1) = 1000.*exp(log(3.*e(1)/(4.*pi))/3.)       ! r: Convert mass [mg] -> radius [um]
    do i = 2, n                                                                            
        e(i) = ax*e(i-1)                            ! Equation (4)
        r(i) = 1000.*exp(log(3.*e(i)/(4.*pi))/3.)
    enddo

    ! initial mass distribution                     !       _
    x0 = xn1/xn0                                    ! x0: L/x^2
    do i = 1, n
        x1   = e(i)
        g(i) = 3.*x1*x1*x0*exp(-x1/xn0)             ! Equation (2), (25)
    enddo

    ! courant numbers
    call courant( n, rq0, dlnr, scal, ax, e, g, r, c, ima )

    ! kernel
    call trkern ( n, isw, g, r, e, ck, ec )

    ! multiply kernel with constant timestep and logarithmic grid distance
    do i = 1, n
        do j = 1, n
            ck(i,j) = ck(i,j)*dt*dlnr               ! ck: Kernel, dlnr: delta y
        enddo
    enddo

    ! collision
    call coad( n, e, g, r, c, ima, ck, ec )

    ! mass balance
    x0 = 0.; x1 = 0.
    do i = 1, n
        x0 = x0 + g(i)*dlnr                 ! Equation (5)
        x1 = max(x1,g(i))
        if ( abs(x1-g(i)) < 1.e-9 ) imax = i
    enddo

    do i = 1, n
        x1   = e(i)
        nr(i) = 1./(3.*x1*x1)*g(i)          ! Equation (2)
    enddo

end subroutine 

subroutine coad( n, e, g, r, c, ima, ck, ec )
    integer,intent(in)    :: n
    real, dimension(n,n), intent(in)    :: ck
    real, dimension(n),   intent(inout) :: e, g, r
    real, dimension(n,n), intent(inout) :: c, ima, ec
    !common /cour/ c(n,n),ima(n,n)
    !common /grid/ g(n),r(n),e(n)
    !common /kern/ ck(n,n),ec(n,n)

    integer :: i0, i1, i, j, k, kp
    real :: gmin, x0, gsi, gsj, gsk, gk, x1, flux 
    data gmin /1.e-60/ 

! collision subroutine, exponential approach
!      parameter (n=400)
!      implicit double precision (a-h,o-z)
!      common /cour/ c(n,n),ima(n,n)
!      common /grid/ g(n),r(n),e(n)
!      common /kern/ ck(n,n),ec(n,n)
!      data gmin /1.d-60/

    ! lower and upper integration limit i0,i1
    do i = 1, n-1
        i0 = i
        if ( g(i) > gmin ) go to 2000
    enddo
2000 continue

    do i = n-1, 1, -1
        i1 = i
        if ( g(i) > gmin ) go to 2010
    enddo
2010 continue

    do i = i0, i1
        do j = i, i1
            k  = ima(i,j)
            kp = k+1
            x0 = ck(i,j)*g(i)*g(j)
            x0 = min(x0,g(i)*e(j))
            if (j /= k) x0 = min(x0,g(j)*e(i))
            gsi  = x0 / e(j)
            gsj  = x0 / e(i)
            gsk  = gsi + gsj
            g(i) = g(i) - gsi
            g(j) = g(j) - gsj
            gk   = g(k) + gsk

            if ( gk > gmin ) then
                x1    = log( g(kp) / gk+1.e-60 )
                flux  = gsk / x1*( exp(0.5*x1)-exp(x1*(0.5-c(i,j))) )
                flux  = min( flux, gk )
                g(k)  = gk - flux
                g(kp) = g(kp) + flux
            endif
        enddo
    enddo

    return
end subroutine coad


subroutine courant ( n, rq0, dlnr, scal, ax, e, g, r, c, ima )
    integer,          intent(in)                 :: n
    real, intent(in)                 :: rq0, dlnr, scal, ax
    real, intent(in), dimension(n)   :: e, g, r
    real, intent(out),dimension(n,n) :: c, ima
    
    integer            :: i, j, k, kk
    real   :: x0  
    !implicit double precision (a-h,o-z)
    !common /const/ rq0,dlnr,scal,ax
    !common /cour/ c(n,n),ima(n,n)
    !common /grid/ g(n),r(n),e(n)
    do i = 1, n
        do j = i, n
            x0 = e(i) + e(j)
            do k = j, n
                if ( (e(k) >= x0) .and. (e(k-1) < x0) ) then
                    if ( c(i,j) < (1.-1.e-08) ) then
                        kk     = k-1
                        c(i,j) = log( x0/e(k-1) ) / ( 3.e0*dlnr )
                    else
                        c(i,j) = 0.
                        kk     = k
                    endif
                    ima(i,j) = min(n-1,kk)
                    go to 2000
                endif
            enddo
        2000 continue
            c(j,i)   = c(i,j)
            ima(j,i) = ima(i,j)
        enddo
    enddo

    return
end subroutine


subroutine trkern( n, isw, g, r, e, ck, ec )
    integer,          intent(in)                   :: isw, n
    real, intent(in),   dimension(n)   :: g,r,e
    real, intent(inout),dimension(n,n) :: ck,ec
    
    integer                          :: i, j, k, jm, im, jp, ip
    real                 :: x0, effi
    real, parameter      :: pi = 3.141592
    real, dimension(n)   :: rr, winf
    real, dimension(n,n) :: cck

    !parameter (n=400)
    !implicit double precision (a-h,o-z)
    !common /grid/ g(n),r(n),e(n)
    !common /kern/ ck(n,n),ec(n,n)
    !common /veloc/ winf(n),rr(n)
    !dimension cck(n,n)
    !data pi/3.141592654/


    ! terminal velocity
    call fallg( n, g, r, e, rr, winf )
    if ( isw == 0 ) then
        ! long kernel
        do j = 1, n
            do i = 1, j
                if( r(j) <= 50. ) then
                    effi = 4.5e-4*r(j)*r(j)*   &
                          ( 1.e0-3.e0/(max(3.e0,(r(i)))+1.e-2) )
                else
                    effi = 1.e0
                endif
            cck(j,i) = pi*(rr(j)+rr(i))*(rr(j)+rr(i))*effi*   &
                       abs(winf(j)-winf(i))
            cck(i,j) = cck(j,i)
            enddo
        enddo
    elseif ( isw == 1 ) then
        ! halll
        call effic( n, g, r, e, ec )
        do j = 1, n
            do i = 1, j
                cck(j,i) = pi*(rr(j)+rr(i))*(rr(j)+rr(i))*ec(j,i)*  &
                           abs(winf(j)-winf(i))
                cck(i,j) = cck(j,i)
            enddo
        enddo
    else
        ! golovin kernel
        do j = 1, n
            do i = 1, j
                cck(j,i) = 1.5*( e(j)+e(i) )
                cck(i,j) = cck(j,i)
            enddo
        enddo
    endif

    ! two-dimensional linear interpolation of kernel
    do i = 1, n
        do j = 1, n
            jm = max0(j-1,1)
            im = max0(i-1,1)
            jp = min0(j+1,n)
            ip = min0(i+1,n)
            ck(i,j) = 0.125*( cck(i,jm) + cck(im,j) + &
                              cck(ip,j) + cck(i,jp) ) + .5*cck(i,j)
            if ( i == j ) ck(i,j) = 0.5*ck(i,j)
        enddo
    enddo

    return
end subroutine


subroutine fallg( n, g, r, e, rr, winf)
! terminal velocity of falling drops
    integer, intent(in)                          :: n
    real, intent(in),dimension(n)    :: g,r,e
    real, intent(out),dimension(n)   :: rr
    real, intent(inout),dimension(n) :: winf

    real, dimension(7)  :: b    
    real, dimension(6)  :: c    
    real, dimension(20) :: rat
    real, dimension(15) :: r0
    integer            :: i,j,k
    real, parameter    :: pi=3.141592654
    real   :: x0, eta, xlamb, rhow, rhoa, grav, cunh, &
                          t0, sigma, stok, stb, phy, py, bond, x,y,xrey  


    !parameter (n=400)
    !implicit double precision (a-h,o-z)
    !common /grid/ g(n),r(n),e(n)
    !common /kern/ ck(n,n),ec(n,n)
    !common /veloc/ winf(n),rr(n)
    !dimension b(7),c(6),rat(20),r0(15),ecoll(15,20)
    data b /-0.318657e1,0.992696,-0.153193e-2,-0.987059e-3, &
           -0.578878e-3,0.855176e-4,-0.327815e-5/
    data c /-0.500015e1,0.523778e1,-0.204914e1,0.475294,-0.542819e-1, &
            0.238449e-2/

    eta   = 1.818e-4
    xlamb = 6.62e-6
    rhow  = 1.
    rhoa  = 1.225e-3
    grav  = 980.665
    cunh  = 1.257*xlamb
    t0    = 273.15
    sigma = 76.1-0.155*(293.15-t0)
    stok  = 2.*grav*(rhow-rhoa)/(9.*eta)
    stb   = 32.*rhoa*(rhow-rhoa)*grav/(3.*eta*eta)
    phy   = sigma*sigma*sigma*rhoa*rhoa/((eta**4)*grav*(rhow-rhoa))
    py    = phy**(1./6.)

    ! rrradius in cm-units
    do j = 1, n
        rr(j) = r(j)*1.e-4
    enddo

    do j = 1, n
        if ( rr(j) <= 1.e-3 ) then
            winf(j) = stok*(rr(j)*rr(j)+cunh*rr(j))
        elseif ( (rr(j) > 1.e-3) .and. (rr(j) <= 5.35e-2) ) then
            x = log(stb*rr(j)*rr(j)*rr(j))
            y = 0.
            do i = 1, 7
                y = y + b(i)*(x**(i-1))
            enddo
            xrey    = ( 1.+cunh/rr(j) ) * exp(y)
            winf(j) = xrey*eta / ( 2.*rhoa*rr(j) )
        elseif ( rr(j) > 5.35e-2 ) then
            bond = grav*(rhow-rhoa)*rr(j)*rr(j)/sigma
            if ( rr(j) > 0.35 ) bond = grav*(rhow-rhoa)*0.35*0.35/sigma
            x = log(16.*bond*py/3.)
            y = 0.
            do i = 1, 6
                y = y + c(i)*(x**(i-1))
            enddo
            xrey    = py*exp(y)
            winf(j) = xrey*eta / ( 2.*rhoa*rr(j) )
            if ( rr(j) > 0.35 ) winf(j) = xrey*eta / ( 2.*rhoa*0.35 )
       endif
    enddo
    return
end subroutine

subroutine effic ( n, g, r, e, ec )
!collision efficiencies of hall kernel
    integer, intent(in)                 :: n
    real, intent(in),    dimension(n)   :: g,r,e
    real, intent(inout), dimension(n,n) :: ec
    real, dimension(21)    :: rat
    real, dimension(15)    :: r0
    real, dimension(15,21) :: ecoll

    ! common /grid/ g(n),r(n),e(n)
    ! common /kern/ ck(n,n),ec(n,n)
    ! dimension rat(21),r0(15),ecoll(15,21)
    data r0 /6.,8.,10.,15.,20.,25.,30.,40.,50., &
             60.,70.,100.,150.,200.,300./
    data rat /0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5, &
              0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0/
    data ecoll / &
     0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001 &
    ,0.001,0.001,0.001,0.001,0.001,0.003,0.003,0.003,0.004,0.005 &
    ,0.005,0.005,0.010,0.100,0.050,0.200,0.500,0.770,0.870,0.970 &
    ,0.007,0.007,0.007,0.008,0.009,0.010,0.010,0.070,0.400,0.430 &
    ,0.580,0.790,0.930,0.960,1.000,0.009,0.009,0.009,0.012,0.015 &
    ,0.010,0.020,0.280,0.600,0.640,0.750,0.910,0.970,0.980,1.000 &
    ,0.014,0.014,0.014,0.015,0.016,0.030,0.060,0.500,0.700,0.770 &
    ,0.840,0.950,0.970,1.000,1.000,0.017,0.017,0.017,0.020,0.022 &
    ,0.060,0.100,0.620,0.780,0.840,0.880,0.950,1.000,1.000,1.000 &
    ,0.030,0.030,0.024,0.022,0.032,0.062,0.200,0.680,0.830,0.870 &
    ,0.900,0.950,1.000,1.000,1.000,0.025,0.025,0.025,0.036,0.043 &
    ,0.130,0.270,0.740,0.860,0.890,0.920,1.000,1.000,1.000,1.000 &
    ,0.027,0.027,0.027,0.040,0.052,0.200,0.400,0.780,0.880,0.900 &
    ,0.940,1.000,1.000,1.000,1.000,0.030,0.030,0.030,0.047,0.064 &
    ,0.250,0.500,0.800,0.900,0.910,0.950,1.000,1.000,1.000,1.000 &
    ,0.040,0.040,0.033,0.037,0.068,0.240,0.550,0.800,0.900,0.910 &
    ,0.950,1.000,1.000,1.000,1.000,0.035,0.035,0.035,0.055,0.079 &
    ,0.290,0.580,0.800,0.900,0.910,0.950,1.000,1.000,1.000,1.000 &
    ,0.037,0.037,0.037,0.062,0.082,0.290,0.590,0.780,0.900,0.910 &
    ,0.950,1.000,1.000,1.000,1.000,0.037,0.037,0.037,0.060,0.080 &
    ,0.290,0.580,0.770,0.890,0.910,0.950,1.000,1.000,1.000,1.000 &
    ,0.037,0.037,0.037,0.041,0.075,0.250,0.540,0.760,0.880,0.920 &
    ,0.950,1.000,1.000,1.000,1.000,0.037,0.037,0.037,0.052,0.067 &
    ,0.250,0.510,0.770,0.880,0.930,0.970,1.000,1.000,1.000,1.000 &
    ,0.037,0.037,0.037,0.047,0.057,0.250,0.490,0.770,0.890,0.950 &
    ,1.000,1.000,1.000,1.000,1.000,0.036,0.036,0.036,0.042,0.048 &
    ,0.230,0.470,0.780,0.920,1.000,1.020,1.020,1.020,1.020,1.020 &
    ,0.040,0.040,0.035,0.033,0.040,0.112,0.450,0.790,1.010,1.030 &
    ,1.040,1.040,1.040,1.040,1.040,0.033,0.033,0.033,0.033,0.033 &
    ,0.119,0.470,0.950,1.300,1.700,2.300,2.300,2.300,2.300,2.300 &
    ,0.027,0.027,0.027,0.027,0.027,0.125,0.520,1.400,2.300,3.000 &
    ,4.000,4.000,4.000,4.000,4.000/

    integer :: i, j, k, ir, kk, iq
    real    :: rq, p, q, ek 

    ! two-dimensional linear interpolation of the collision efficiency
    do j = 1, n
        do i = 1, j
            do k = 2, 15
                if ( (r(j) <= r0(k)) .and. (r(j) >= r0(k-1)) ) then
                    ir = k
                elseif ( r(j) > r0(15) ) then
                    ir = 16
                elseif ( r(j) < r0(1) ) then
                    ir = 1
                endif
            enddo

            rq = r(i)/r(j)

            do kk = 2, 21
                if ( (rq <= rat(kk)) .and. (rq > rat(kk-1)) ) iq = kk
            enddo

            if ( ir < 16 ) then
                if ( ir >= 2 ) then
                    p = ( r(j)-r0(ir-1) ) / ( r0(ir)-r0(ir-1) )
                    q = ( rq-rat(iq-1) ) / ( rat(iq)-rat(iq-1) )
                    ec(j,i) = (1.-p)*(1.-q)*ecoll(ir-1,iq-1)+ &
                              p*(1.-q)*ecoll(ir,iq-1)+        &
                              q*(1.-p)*ecoll(ir-1,iq)+        &
                              p*q*ecoll(ir,iq)
                else
                    q = ( rq-rat(iq-1) ) / ( rat(iq)-rat(iq-1) )
                    ec(j,i) = (1.-q)*ecoll(1,iq-1)+q*ecoll(1,iq)
                endif
            else
                q       = ( rq-rat(iq-1) ) / ( rat(iq)-rat(iq-1) )
                ek      = (1.-q)*ecoll(15,iq-1)+q*ecoll(15,iq)
                ec(j,i) = min(ek,1.0)
            endif

            ec(i,j) = ec(j,i)

            if ( ec(i,j) < 1.e-20 ) stop 99
        enddo
    enddo

    return
end subroutine

end module collision_mod
