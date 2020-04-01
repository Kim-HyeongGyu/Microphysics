module advection_mod
contains
    subroutine compute_advection(w_full, initC, sfc_C, nt, nz, &
                                 z_full, z_half, diff_method, C)
    ! Vertical advection = w * (dC/dz)
    ! diff_method: 1. Foward/Centerd finite difference
    !              2. Finite volume method
    implicit none
    integer, intent(in)    :: nt
    character(len=*), intent(in) :: diff_method
    integer :: n, k, nz, dt = 10      ! CFL condition (mu=w*dt/dz)
    real    :: dC_dz, dz, zbottom = 0.
    real    :: sfc_C    ! for surface effect (radiation, etc ...)
    real    :: flux_minus, flux_plus   ! for finite volume method
    real, dimension(nz)    :: w_full, initC, z_full
    real, dimension(nz+1)  :: w_half, z_half
    real, dimension(nt,nz) :: C

    ! Initial value
    C(1,:) = initC

    ! Make stagged grid for advection
    w_half(2:nz) = ( w_full(2:) + w_full(1:nz-1) ) / 2.

    ! homogeneous Dirichlet Boundary Condition
    w_half(1) = 0.; w_half(nz+1) = 0.  ! Top/Bottom w=0  [m s-1]

    select case (diff_method)
        ! Here, we use Lorenz configuration
        ! See Figure 1 in Holdaway et al., (2012) 
        ! https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1002/qj.2016
        case ("finite_difference")
            !  C(n+1,k) - C(n,k)           C(n,k+1) - C(n,k-1)
            ! ------------------- = -w(k) ---------------------
            !          dt                  z(n,k+1) - z(n,k-1)
            do n = 1, nt-1
                do k = 2, nz-1
                    dC_dz = ( C(n,k+1) - C(n,k-1) ) / &
                            ( z_full(k+1) - z_full(k-1) )
                    C(n+1,k) = C(n,k) - w_full(k)*dt*dC_dz
                end do

                ! Foward/Backward difference at boundary
                C(n+1, 1) = C(n,1) - w_half(2)*dt* &
                          ( C(n,2) - C(n,1) ) /  &
                          ( z_full(2) - z_full(1) )
                C(n+1,nz) = C(n,nz) - w_half(nz)*dt* &
                          ( C(n,nz) - C(n,nz-1) ) /  &
                          ( z_full(nz) - z_full(nz-1) )
            end do
        case ("finite_volume")
            do n = 1, nt-1
                ! TODO: Need to discuss BC
                ! Boundary condition
                k = 1
                flux_minus = sfc_C *w_half(k)*dt    ! equal zero
                flux_plus  = C(n,k)*w_half(k+1)*dt
                dz         = z_half(k+1) - z_half(k)
                C(n+1,k)   = C(n,k)+(flux_minus-flux_plus)/dz

                do k = 2, nz
                    flux_minus = C(n,k-1)*w_half(k)*dt
                    flux_plus  = C(n,k  )*w_half(k+1)*dt
                    dz         = z_half(k+1) - z_half(k)
                    C(n+1,k)   = C(n,k)+(flux_minus-flux_plus)/dz
                end do
            end do
        case default
            print*, "Not setup diff_method option. &
                     please check input.nml"
            stop
    end select

    end subroutine compute_advection
end module advection_mod
