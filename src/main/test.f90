    program test
    
    implicit none

    real :: L, Rv, Ka, Temp, Fk, rho, Dv, Pres, es, Fd

    Temp = 273.15   ! [K]
    Pres = 700      ! [hPa]

    L  = 2500297.8  ! heat of vaporization                          [J kg-1] 
    Rv = 467        ! vapor gas constant                            [J kg-1 K-1]
    rho= 1000.      ! water density                                 [kg m-3]

    ! Refer to Rogers & Yau (1996), 103p - Table 7.1 (T=273K)
    Ka = 2.40e-2    ! coefficient of thermal conductivity of air    [J m-1 s-1 K-1]
    Dv = 2.21e-5    ! molecular diffusion coefficient               [m2 s-1]


    Fk = ((L/(Rv*Temp))-1.)*((L*rho)/(Ka*Temp))
    print*, 'Fk = ',Fk

    es = 6.112 * exp(( 17.67*(Temp-273.15) )/( (Temp-273.15)+243.5 )) 
    ! To calculate Fd, need to convert the units of 'es'. :: [hPa] > [J m-3]
    Fd = (rho*Rv*Temp) / ((Dv*(1000./Pres))*(es*100.))
    print*, Dv*(1000./Pres), es
    print*, 'Fd = ',Fd
   
    
    
    endprogram
