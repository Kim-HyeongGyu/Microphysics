module file_io_mod
use            global_mod
use            read_nc_mod
use            write_nc_mod
implicit none
    
contains

    subroutine read_namelist()  ! {{{
        ! Set variables
        namelist /main_nml/     t_final
        namelist /data_nml/     vert_var, temp_var
        namelist /dynamics_nml/ num_levels, top_of_atmosphere,  &
                                vertical_grid, vertical_advect, &
                                CFL_condition, status_case
        namelist /physics_nml/  drop_var, rmin, rmax, rratio,   &
                                nbin, Nc, qc, dist_type

        open  (unit = 8, file = 'input.nml', delim = 'apostrophe')
        read  (unit = 8, nml  = main_nml) 
        read  (unit = 8, nml  = data_nml) 
        read  (unit = 8, nml  = dynamics_nml) 
        read  (unit = 8, nml  = physics_nml) 
        close (unit = 8)

        ! Overide default values for optional arguments
        nz      = num_levels
        ztop    = top_of_atmosphere   ! [m]
    end subroutine read_namelist    ! }}}



    subroutine read_data_init(nlev,lev,temp_in,qv_in,w_in) ! {{{

    integer,intent(out) :: nlev
    real, dimension(:),allocatable,intent(out) :: lev, temp_in, qv_in, w_in
    character(len=5), parameter :: vname1 = "t"
    character(len=5), parameter :: vname2 = "q"
    character(len=5), parameter :: vname3 = "w"
    character(len=300) :: INAME1, INAME2, INAME3, IPATH
    character(len=6) :: ld_name = "level", l_name= "lev" ! lev dimension_name, lev_var_name

     IPATH = './input'

     WRITE(INAME1,'(4A)') trim(IPATH),'/ERA_Interim_',trim(vname1),'_AUG.nc'
     WRITE(INAME2,'(4A)') trim(IPATH),'/ERA_Interim_',trim(vname2),'_AUG.nc'
     WRITE(INAME3,'(4A)') trim(IPATH),'/ERA_Interim_',trim(vname3),'_AUG.nc'

     call griddims(INAME1,ld_name,nlev)

      allocate(lev(nlev), temp_in(nlev), qv_in(nlev),w_in(nlev))
     call read_nc_data(INAME1,vname1,l_name,nlev,lev,temp_in)
     call read_nc_data(INAME2,vname2,l_name,nlev,lev,qv_in)
     call read_nc_data(INAME3,vname3,l_name,nlev,lev,w_in)

        ! Vertical wind [m s-1]
      w_in      = sin( (/ (I, I = 1,nlev*2,2) /) / 10. )
!        w_in    = w_in*0. + 2.

!      ========== ideal status ==============
       if (status_case == "ideal") then
        ! Temperature   [K]
        temp_in = (/ (I+273, I = nlev,1,-1) /)  ! lapse rate 1K/km

        ! Mixing ratio  [kg kg-1]
        qv_in(:) = w_in*0.
        qv_in(5) = 0.5
        end if     


    end subroutine read_data_init   ! }}}


    subroutine write_data() ! {{{

    character(len=300) :: ONAME1, ONAME2, OPATH
    character(len=5), parameter :: vname1 = "T"
    character(len=5), parameter :: vname2 = "q"

     OPATH = './output'
     WRITE(ONAME1,'(10A)') trim(OPATH),'/',trim(vname1),'_',trim(status_case),'_',trim(vertical_grid),'_',trim(vertical_advect),'_out.nc'
     WRITE(ONAME2,'(10A)') trim(OPATH),'/',trim(vname2),'_',trim(status_case),'_',trim(vertical_grid),'_',trim(vertical_advect),'_out.nc'
     call write_nc_data(ONAME1,vname1,nt,nz,z_full,T)
     call write_nc_data(ONAME2,vname2,nt,nz,z_full,q)


        !open(unit = 30, file = "Tout.txt")
        ! open(unit = 30, file = "output.dat", form='unformatted', &
        !      status = "unknown", access='direct', recl=4*nz) 
        !do n = 1, nt
        !    write(30,*) T(:,n)
        !end do
        !close(30)

        !open(unit = 30, file = "qout.txt")
        ! open(unit = 30, file = "output.dat", form='unformatted', &
        !      status = "unknown", access='direct', recl=4*nz) 
        !do n = 1, nt
        !    write(30,*) q(:,n)
        !end do
        !close(30)
    end subroutine write_data   ! }}}

end module file_io_mod
