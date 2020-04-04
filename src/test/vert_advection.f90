module advection_mod
contains
    subroutine compute_advection(w_full, initC, sfc_C, nt, nz, &
                                 dz, diff_method, C)
    ! Vertical advection = w * (dC/dz)
    ! diff_method: 1. Foward/Centerd finite difference
    !              2. Finite volume method
    implicit none
    integer, intent(in)    :: nt
    character(len=*), intent(in) :: diff_method
    integer :: n, k, nz, dt = 10      ! CFL condition (mu=w*dt/dz)
    real    :: dC_dz, zbottom = 0.
    real    :: sfc_C    ! for surface effect (radiation, etc ...)
    real    :: flux_minus, flux_plus   ! for finite volume method
    ! real, dimension(nz)    :: w_full, initC, z_full, dz 
    ! real, dimension(nz+1)  :: w_half, z_half, C_half
    real, dimension(nz)    :: w_full, initC, dz 
    real, dimension(nz+1)  :: w_half, C_half, flux
    real, dimension(nz,nt) :: C

    ! Initial value
    C(:,1) = initC

    ! Make stagged grid for advection
    w_half(2:nz) = ( w_full(2:) + w_full(1:nz-1) ) / 2.

    ! homogeneous Dirichlet Boundary Condition
    w_half(1) = 0.; w_half(nz+1) = 0.  ! Top/Bottom w=0  [m s-1]

    select case (diff_method)
        ! Here, we use Lorenz configuration
        ! See Figure 1 in Holdaway et al., (2012) 
        ! https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1002/qj.2016

        ! 1. Finite difference method {{{
        case ("finite_difference")
            !  C(n+1,k) - C(n,k)           C(n,k+1) - C(n,k-1)
            ! ------------------- = -w(k) ---------------------
            !          dt                  z(n,k+1) - z(n,k-1)
            do n = 1, nt-1
                do k = 2, nz-1
                    dC_dz = ( C(k+1,n) - C(k-1,n) ) / &
                            ( dz(k) + (dz(k+1) + dz(k-1))/2. )
                            ! ( z_full(k+1) - z_full(k-1) )
                    C(k,n+1) = C(k,n) - w_full(k)*dt*dC_dz
                end do

                ! Foward/Backward difference at boundary
                C( 1,n+1) = C(1,n) - w_half(2)*dt* &
                          ( C(2,n) - C(1,n) ) /  &
                          ( (dz(2)+dz(1))/2. )
                          ! ( z_full(2) - z_full(1) )
                C(nz,n+1) = C(nz,n) - w_half(nz)*dt* &
                          ( C(nz,n) - C(nz-1,n) ) /  &
                          ( (dz(nz)+dz(nz-1))/2. )
                          ! ( z_full(nz) - z_full(nz-1) )
            end do !}}}

        ! 2) Finite Volume method (FVM) {{{
        case ("finite_volume") 
            do n = 1, nt-1
                do k = 2, nz
                    C_half(k) = ( dz(k-1)*C(k-1,n) + dz(k)*C(k,n) ) &
                              / ( dz(k-1) + dz(k) )
                end do

                ! Boundary values
                C_half(1)    = C( 1,n) - ( C_half( 2)-C( 1,n) )
                C_half(nz+1) = C(nz,n) - ( C_half(nz)-C(nz,n) )

                flux = w_half*C_half*dt
                do k = 1, nz
                    C(k,n+1) = C(k,n) + ( flux(k)-flux(k+1) ) / dz(k)
                end do
            end do !}}}

        case default
            print*, "Not setup diff_method option. &
                     please check input.nml"
            stop
    end select

    end subroutine compute_advection
end module advection_mod
