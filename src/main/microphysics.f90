    program bin_microphysics

    implicit none

    character(len=*), intent(in) :: dist_type
    integer :: nz, nbin = 40
    real*8, dimension(nbin) :: N
    real, dimension(nbin)   :: r , m, dr
    real, dimension(nbin+1) :: rb, mb
    real :: ratio=1.3   ! common ratio
    real :: rmin=1E-6   ! [m]
    real :: PI=3.141592
    real :: N0 = 1.
    real :: mu, lamda
    real*8 :: ave, std
    integer :: i
    
    rb(1)=rmin
    do i = 1, nbin
        rb(i+1) = rb(i)*ratio
    enddo
!    do i = 1, nbin+1                   ! TODO :: NEED MODIFY ! 
!        mb(i) = (4./3.)*PI*rb(i)**3    !      :: m(i) = N(i)*(density of water)*(4/3)*(PI)*(r(i)^3)
!    enddo
    do i = 1, nbin
        r(i)=(rb(i)+rb(i+1))/2.
!        m(i)=(mb(i)+mb(i+1))/2.
    enddo

    select case (dist_type)
        case ("log_normal_dist")
            ave = 1. ; std = 0.
            !ave = 1. ; std = 4.
            do i = 1, nbin
                ave = ave*r(i)
            enddo
            ave = log(ave**(1./nbin))
            do i = 1, nbin
                std = std + (log(r(i)-ave)**2)
            enddo
                std = sqrt(std/nbin)
            do i = 1, nbin
                !N(i) = (N0/(sqrt(2.*PI)*log(std)*r(i)))*(exp(-(((log(r(i))-log(ave))**2)/(2.*((log(std))**2)))))
                N(i) = (N0/(sqrt(2.*PI)*std))*(exp(-((log(r(i))-ave)**2)/(2.*(std)**2)))
                print*, i, N(i)
            enddo
   
            OPEN(10, FILE='log-normal_dist.gdat', FORM='UNFORMATTED', ACCESS='DIRECT', &
                    STATUS='UNKNOWN', RECL=8)
            do i =1, nbin
                WRITE(10, REC=i) N(i)
            enddo
        case ("gamma_dist")
    
            !mu = 0.1 ; lamda = 1.
            !do i = 1, nbin
            !N(i) = N0*(r(i)**mu)*(exp(-(lamda*r(i))))
            !print*, i, N(i)
            !enddo
    
            !OPEN(10, FILE='gamma_dist.gdat', FORM='UNFORMATTED', ACCESS='DIRECT', &
            !         STATUS='UNKNOWN', RECL=4)
            !do i =1, nbin
            !   WRITE(10, REC=i) N(i)
            !enddo
    endprogram
