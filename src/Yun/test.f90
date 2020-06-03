    program test
    
    implicit none

    integer, dimension(3,4) :: aa = reshape([11,12,13,14,15,16,17,18,19,20,21,22],[3,4])
    integer                 :: i, j

    do i = 1, 3
    do j = 1, 4
        print*, i, j, aa(i,j)
    enddo 
    enddo 
    
    print*, " "
    print*, aa
    
    
    endprogram
