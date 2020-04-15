      program read_nc
      use read_nc_mod


      implicit none
      
      real, dimension(:),allocatable :: lev, var
      integer :: nlev
      character(len=5), parameter :: vname1 = "q"
      character(len=5), parameter :: vname2 = "t"
      character(len=5), parameter :: vname3 = "w"
      character(len=300) :: FNAME1, FNAME2, FNAME3
      character(len=6) :: ld_name = "level", l_name= "lev" ! lev dimension_name, lev_var_name

      WRITE(FNAME1,'(3A)') './ERA_Interim_',trim(vname1),'_AUG.nc'
      WRITE(FNAME2,'(3A)') './ERA_Interim_',trim(vname2),'_AUG.nc'
      WRITE(FNAME3,'(3A)') './ERA_Interim_',trim(vname3),'_AUG.nc'
      print*, FNAME1

      call griddims(FNAME1,ld_name,nlev)

      allocate(lev(nlev))
      allocate(var(nlev))
      print*, nlev
      call read_nc_data(FNAME1,vname1,l_name,nlev,lev,var)
      print*,lev(:)
      print*, var(:) 


      end program read_nc
