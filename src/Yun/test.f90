    program test
    
    implicit none
    
    integer :: i 
    integer :: nbin = 40
    real*8, dimension(nbin)  :: N

    open(10, file='log-normal_dist.gdat', form='unformatted', access='direct', &
             status='old', recl=8)
    do i = 1, nbin
        READ(10,rec=i) N(i)
        print*, N(i)
    enddo

    endprogram
