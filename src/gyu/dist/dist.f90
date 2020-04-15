program dist
use microphysics

    implicit none
    integer :: i
    real :: NC = 1.e2   ! 10~1000
    real :: qc = 2.     ! 0.5~2
    real, dimension(40) :: radius, Nr
    ! double precision, dimension(100) :: Nr
    character(len=30) :: dist_type
    integer :: nbin=40
    real :: rratio=1.26  ! common ratio
    real :: rmin=1.e-6   ! [m] = 1 [um]
    
    dist_type = "gamma"

    call make_bins(rmin, rratio, nbin, radius)
    print*, radius
    call conc_dist(radius, Nc, qc, dist_type, Nr)
    print*, dist_type

    open(unit = 30, file = "Nr.txt")
    write(30,*) radius
    write(30,*) Nr
    close(30)
    
end program dist

