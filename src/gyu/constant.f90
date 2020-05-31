module constant_mod
! Constant module
! all variable units are MKS based.
implicit none
    
    ! Constant parameter
    real, parameter :: P0  = 1000.00        ! [hPa] Surface pressure
    real, parameter :: PI  = 3.141592
    real, parameter :: R   = 287.           ! [J kg-1 K-1] Gas constant
    real, parameter :: Cp  = 1005.          ! [J kg-1 K-1] Specific heat 
                                            !             at constant pressure
    real, parameter :: g   = 9.8            ! [m s-2] Gravity 

    real, parameter :: rho_liquid = 1000.   ! [kg m-3] Water density
    real, parameter :: rho_air    = 1.225   ! [kg m-3] Air density

    real, parameter :: L  = 2500297.8       ! [J kg-1] Heat of vaporization 
    real, parameter :: Rv = 467             ! [J kg-1 K-1] Vapor gas constant   

    ! Refer to Rogers & Yau (1996), 103p - Table 7.1 (T=273K)
    real, parameter :: Ka = 2.40e-2         ! [J m-1 s-1 K-1] Coefficient of thermal 
                                            !                conductivity of air
    real, parameter :: Dv = 2.21e-5         ! [m2 s-1] Molecular diffusion coefficient
    real, parameter :: mu = 1.717e-5        ! [kg m-1 s-1] Dynamic viscosity (at 273 [K])
    real, parameter :: sigma = 7.5e-2       ! [N m-1] Surface tension (Yau (1996) - 85p)

end module constant_mod
