       program micro_test


       implicit none


       integer :: i
       real, parameter :: pi = 3.141592
       integer, parameter :: nm = 40
       real, parameter :: ratio = 1.3, n0=1.

       real, dimension(nm) :: r, nor
       real, dimension(nm+1) :: rb
       real :: mm, std,ave,mu,sigma


!       mass = 4/3*pi*r**3
        

       rb(1) = 1E-6
     
       std = 2.
       ave = 1.


!       m = (4./3.)*pi*(r(1))**3

       do i = 2, nm+1

!         rb(i) = rb(1)*(ratio**3)**(i-1)

          rb(i) = rb(1)*((ratio)**(i-1))
          r(i-1) = (rb(i-1)+rb(i))/2.
!       write(*,'(2E12.4)') rb(i-1), r(i-1)
!         print*,rb(i)
       end do

 
!       do i = 1, nm
!        print*,i,r(i) ,log(r(i))

!       end do

!       mu = log(ave**2/sqrt(std+ave**2))
!       sigma = sqrt(log(std/(ave**2+1)))


!========================================================
       do i = 1, nm
!       nor(i)=n0/((sqrt(2*pi))*sigma*r(i))*(exp(-(log(r(i)-mu))**)/2*sigma**2)
       nor(i)=(n0/((sqrt(2*pi))*log(std)*r(i)))*(exp(-((log(r(i))-log(ave))**2)/(2*log(std)**2)))
       print*, nor(i)
       end do



       open(90,FILE='./nor_dis.gdat',STATUS='UNKNOWN',FORM='UNFORMATTED', &
               ACCESS='DIRECT', RECL=4)
       do i = 1, nm
       write(90,rec=i) nor(i)
       end do





       end program 
