module namelist_mod
! Variable from namelist
integer :: nz             = 20             ! [#] Num of levels
integer :: nt             = 10000          ! [#] iteration number of time
integer :: nbin           = 100            ! [#] Number of bins  
integer :: header_data    = 0              ! header line from input
integer :: drop_var       = 1              ! 1: rmin & rmax, 2: rmin & rratio
integer :: grid_type      = 1              ! 1: Constant grid
!          grid_type      = 2              ! 2: stretching grid

!          dyn_adv_scheme = 1              ! 1: Finite difference
integer :: dyn_adv_scheme = 2              ! 2: Finite Volume
!          dyn_adv_scheme = 3              ! 3: PPM (Colella and Woodward, 1984)
!          dyn_adv_scheme = 4              ! 4: PPM (Lin, 2003)

integer :: phy_adv_scheme = 1              ! 1: Reassign
!          phy_adv_scheme = 2              ! 2: Finite Volume
!          phy_adv_scheme = 3              ! 3: PPM (Colella and Woodward, 1984)

integer :: dist_type      = 1              ! 1: Log-normal distribution
!          dist_type      = 2              ! 2: Gamma distribution

real    :: dt      = 1.                    ! [s] time step  
real    :: zbottom = 0.                    ! [m] model bottom
real    :: ztop    = 1000.                 ! [m] model top
real    :: rmin    = 1.e-6                 ! [m] minimum radius of bins 
real    :: rmax    = 1.e-2                 ! [m] maximum radius of bins 
real    :: rratio  = 1.05                  ! [#] radius ratio of bins   
real    :: Nc      = 1.e8                  ! [# m-3] number of droplets     
real    :: qc      = 2.e-3                 ! [kg kg-1] cloud mass mixing ratio
real    :: w_speed = 0.5                   ! [m s-1] vertical wind speed
real    :: w_zero_time = -1                ! [s] Time at zero vertical wind

logical :: ventilation_effect = .false.    ! ventilation effect
logical :: collision_effect   = .false.    ! collision effect
logical :: surface_data       = .false.
logical :: w_from_data        = .false.
! Note! if surface_data = .true. 
!       -> 1st line must be pressure, temperature, mixing ratio, (w)
! Note! if w_from_data = .true. -> w_speed will be updated from data

character(len=20) :: file_name
character(len=20) :: pres_units = "hPa"
character(len=20) :: temp_units = "K"
character(len=20) :: qv_units   = "kg/kg"
character(len=10) :: vert_var
character(len=10) :: temp_var

contains

    subroutine read_namelist()  ! {{{
        namelist /main_nml/     dt,                 &
                                nt
        namelist /data_nml/     file_name,          &
                                header_data,        &
                                surface_data,       &
                                w_from_data,        &
                                pres_units,         &
                                temp_units,         &
                                w_speed,            &
                                qv_units,           & 
                                vert_var,           &
                                temp_var
        namelist /dynamics_nml/ nz,                 &
                                ztop,               &
                                grid_type,          &
                                w_zero_time,        &
                                dyn_adv_scheme
        namelist /physics_nml/  drop_var,           &
                                rmin,               &
                                rmax,               &
                                rratio,             &
                                nbin,               &
                                Nc,                 &
                                qc,                 &
                                dist_type,          &
                                phy_adv_scheme,     &
                                ventilation_effect, &
                                collision_effect

        open  (unit = 8, file = 'input.nml', delim = 'apostrophe')
        read  (unit = 8, nml  = main_nml) 
        read  (unit = 8, nml  = data_nml) 
        read  (unit = 8, nml  = dynamics_nml) 
        read  (unit = 8, nml  = physics_nml) 
        close (unit = 8)

    end subroutine read_namelist    ! }}}
    
    subroutine write_model_setup_info() ! {{{
        character(len=16) :: vartochar

        open (10, file='nml_info.nml') 

        write(10,*), "&main_nml"
        write(10,*), " dt                 = ", dt
        write(10,*), " nt                 = ", nt
        write(10,*), " "

        write(10,*), "&data_nml"
        write(10,*), " file_name          = ", file_name
        write(10,*), " header_data        = ", header_data
        write(10,*), " surface_data       = ", surface_data
        write(10,*), " w_from_data        = ", w_from_data
        write(10,*), " pres_units         = ", pres_units 
        write(10,*), " temp_units         = ", temp_units 
        write(10,*), " w_speed            = ", w_speed
        write(10,*), " qv_units           = ", qv_units    
        write(10,*), " vert_var           = ", vert_var   
        write(10,*), " temp_var           = ", temp_var
        write(10,*), " "
        
        write(10,*) "&dynamics_nml"
        write(10,*) " nz                 = ", nz         
        write(10,*) " ztop               = ", ztop       
        write(10,*) " grid_type          = ", grid_type  
        write(10,*) " dyn_adv_scheme     = ", dyn_adv_scheme
        write(10,*) " "

        write(10,*) "&physics_nml"
        write(10,*) " drop_var           = ", drop_var   
        write(10,*) " rmin               = ", rmin       
        write(10,*) " rmax               = ", rmax       
        write(10,*) " rratio             = ", rratio     
        write(10,*) " nbin               = ", nbin      
        write(10,*) " Nc                 = ", Nc         
        write(10,*) " qc                 = ", qc         
        write(10,*) " dist_type          = ", dist_type
        write(10,*) " phy_adv_scheme     = ", phy_adv_scheme
        write(10,*) " ventilation_effect = ", ventilation_effect
        
        close(10)
    end subroutine write_model_setup_info   ! }}}
        
end module namelist_mod
