module file_io_mod
use            global_mod
use            read_nc_mod
implicit none
    
contains

    subroutine read_namelist()  ! {{{
        ! Set variables
        namelist /main_nml/     t_final
        namelist /dynamics_nml/ num_levels, top_of_atmosphere,  &
                                vertical_grid, vertical_advect, &
                                CFL_condition
        ! namelist /physics_nml/  rmin, ...

        open  (unit = 8, file = 'input.nml', delim = 'apostrophe')
        read  (unit = 8, nml  = main_nml) 
        read  (unit = 8, nml  = dynamics_nml) 
        close (unit = 8)

        ! Overide default values for optional arguments
        nz      = num_levels
        ztop    = top_of_atmosphere   ! [m]
    end subroutine read_namelist    ! }}}



    subroutine read_data_init(nlev,lev,temp_in,qv_in,w) ! {{{

    integer,intent(out) :: nlev
    real, dimension(:),allocatable,intent(out) :: lev, temp_in, qv_in, w
    real, dimension(:),allocatable :: w_in
    character(len=5), parameter :: vname1 = "t"
    character(len=5), parameter :: vname2 = "q"
    character(len=5), parameter :: vname3 = "w"
    character(len=300) :: FNAME1, FNAME2, FNAME3, path
    character(len=6) :: ld_name = "level", l_name= "lev" ! lev dimension_name, lev_var_name

    path = '/home/ojs9294/class/mirco_2020/Microphysics/exp/input'

     WRITE(FNAME1,'(4A)') trim(path),'/ERA_Interim_',trim(vname1),'_AUG.nc'
     WRITE(FNAME2,'(4A)') trim(path),'/ERA_Interim_',trim(vname2),'_AUG.nc'
     WRITE(FNAME3,'(4A)') trim(path),'/ERA_Interim_',trim(vname3),'_AUG.nc'

     call griddims(FNAME1,ld_name,nlev)

      allocate(lev(nlev), temp_in(nlev), qv_in(nlev),w_in(nlev))
     call read_nc_data(FNAME1,vname1,l_name,nlev,lev,temp_in)
     call read_nc_data(FNAME2,vname2,l_name,nlev,lev,qv_in)
     call read_nc_data(FNAME3,vname3,l_name,nlev,lev,w_in)
      print*,lev(:), temp_in(:), qv_in(:), w_in(:)


!      ========== ideal status ==============
!        allocate(Tinit(nz), w(nz), qinit(nz))
        allocate(w(nlev))
        ! Temperature   [K]
!        Tinit(:) = (/ (I+273, I = nz,1,-1) /)  ! lapse rate 1K/km

        ! Vertical wind [m s-1]
!        w      = sin( (/ (I, I = 1,nlev*2,2) /) / 10. )
        w      = w*0. + 2.

        ! Mixing ratio  [kg kg-1]
!       qinit(:) = w*0.
!       qinit(5) = 10.
    end subroutine read_data_init   ! }}}


    subroutine write_data() ! {{{
        open(unit = 30, file = "Tout.txt")
        ! open(unit = 30, file = "output.dat", form='unformatted', &
        !      status = "unknown", access='direct', recl=4*nz) 
        do n = 1, nt
            write(30,*) T(:,n)
        end do
        close(30)

        open(unit = 30, file = "qout.txt")
        ! open(unit = 30, file = "output.dat", form='unformatted', &
        !      status = "unknown", access='direct', recl=4*nz) 
        do n = 1, nt
            write(30,*) q(:,n)
        end do
        close(30)
    end subroutine write_data   ! }}}

end module file_io_mod
