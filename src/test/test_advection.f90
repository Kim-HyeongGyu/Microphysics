module advection_mod
contains
    subroutine compute_advection(w, initC, nt, diff_method, C)
    ! Vertical advection = w * (dC/dz)
    ! diff_method: 1. Foward/Centerd finite difference
    !              2. Finite volume method
    implicit none
    real, intent(in)    :: w, initC, dt
    real, intent(out), dimension(nt,nz) :: C

    ! Initial value
    C(1,nz) = initC

    ! homogeneous Dirichlet Boundary Condition
    w(1) = 0.; w(nz) = 0.

    select case (diff_method)
        case ("finite_difference")
            !  C(n+1,k) - C(n,k)           C(n,k+1) - C(n,k-1)
            ! ------------------- = -w(k) ---------------------
            !          dt                  z(n,k+1) - z(n,k-1)
            do n = 1, nt
                do k = 2, nz-1
                    dC_dz = ( C(n,k+1) - C(n,k-1) ) / &
                            ( z_full(k+1) - z_full(k-1) )
                    C(n+1,k) = C(n,k) - w(k)*dt*dC_dz
                end do
                
                ! Foward/Backward difference at boundary
                C(n+1, 1) = C(n, 1) - w( 1)*dt*( C(n,2) - C(n,1) ) &
                                       / ( z_full(2) - z_full(1) )
                C(n+1,nz) = C(n,nz) - w(nz)*dt*( C(n,nz) - C(n,nz-1) ) &
                                       / ( z_full(nz) - z_full(nz-1) )
            end do
        case ("finite_volume")
        ! Here, we use Lorenz configuration
        ! See Figure 1 in Holdaway et al., (2012) 
        ! https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1002/qj.2016
            print*, "TODO! need to work..."
            stop
            ! flux_minus = C(j-1)*w(j-0.5)*dt
            ! flux_plus  = C(j  )*w(j+0.5)*dt
            ! C(j,n+1)   = C(j,n)+(flux_minus-flux_plus)
        case default
            print*, "Not setup diff_method option. &
                     please check input.nml"
    end select

    end subroutine compute_vert_coord
end module advection_mod
