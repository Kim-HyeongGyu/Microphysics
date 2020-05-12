        program dist_test


        real, dimension(6) :: m, n
        integer :: nn, i


        m = (/1.,2.,4.,8.,16.,32./)
        m = (/1.1,2.2,4.4,8.8,17.6,35.2/)

        n = (/0.,3.,30.,5.,1.,1./)
        nn = 40


        do i = 1, 6
        m = m*1.1

        print*, m
        end do













        end program dist_test


