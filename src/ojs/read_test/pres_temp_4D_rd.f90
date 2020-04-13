! This is part of the netCDF package.
! Copyright 2006 University Corporation for Atmospheric Research/Unidata.
! See COPYRIGHT file for conditions of use.

program pres_temp_4D_rd
  use netcdf
  implicit none
  
  ! This is the name of the data file we will read.
  character (len = *), parameter :: FILE_NAME = "./data/ERA-Interim_T_2018.nc"
  integer :: ncid

  ! We are reading 4D data

  integer, parameter :: NDIMS = 4, NRECS = 12
  integer, parameter :: NLVLS = 37, NLATS = 181, NLONS = 360
  character (len = *), parameter :: LVL_NAME = "level"
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  character (len = *), parameter :: REC_NAME = "time"
  integer :: lvl_dimid, lon_dimid, lat_dimid, rec_dimid

  ! The start and count arrays will tell the netCDF library where to
  ! read our data.
  integer :: start(NDIMS), count(NDIMS)

  ! In addition to the latitude and longitude dimensions, we will also
  ! create latitude and longitude variables which will hold the actual
  ! latitudes and longitudes. Since they hold data about the
  ! coordinate system, the netCDF term for these is: "coordinate
  ! variables."
  real :: lats(NLATS), lons(NLONS)
  integer :: lon_varid, lat_varid

  ! We will read surface temperature and pressure fields. In netCDF
  ! terminology these are called "variables."
  character (len = *), parameter :: PRES_NAME="pressure"
  character (len = *), parameter :: TEMP_NAME="t"
  integer :: pres_varid, temp_varid
  integer :: dimids(NDIMS)

  ! We recommend that each variable carry a "units" attribute.
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: PRES_UNITS = "hPa", TEMP_UNITS = "celsius"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"

  ! Program variables to hold the data we will read in. We will only
  ! need enough space to hold one timestep of data; one record.
  real :: pres_in(NLONS, NLATS, NLVLS)
  real :: temp_in(NLONS, NLATS, NLVLS)
  real, parameter :: SAMPLE_PRESSURE = 900.0
  real, parameter :: SAMPLE_TEMP = 9.0

  ! Use these to calculate the values we expect to find.
  real, parameter :: START_LAT = 90.0, START_LON = 0.0

  ! Loop indices
  integer :: lvl, lat, lon, rec, i

  ! To check the units attributes.
  character*80 pres_units_in, temp_units_in
  character*80 lat_units_in, lon_units_in

  ! Open the file. 
  call check( nf90_open(FILE_NAME, nf90_nowrite, ncid) )
  print*, "========= Open the file.=================== "

  ! Get the varids of the latitude and longitude coordinate variables.
  call check( nf90_inq_varid(ncid, LAT_NAME, lat_varid) )
  call check( nf90_inq_varid(ncid, LON_NAME, lon_varid) )
  print*,"========= Get the varids of the latitude and longitude coordinate variables."

  ! Read the latitude and longitude data.
  call check( nf90_get_var(ncid, lat_varid, lats) )
  call check( nf90_get_var(ncid, lon_varid, lons) )
  print*, "=============== Read the latitude and longitude data.======================="
  
  ! Check to make sure we got what we expected.
  do lat = 1, NLATS
     if (lats(lat) /= START_LAT - (lat - 1) * 1.0) stop 2
  end do
  do lon = 1, NLONS
     if (lons(lon) /= START_LON + (lon - 1) * 1.0) stop 2
  end do
  print*,"============ Check to make sure we got what we expected.==================="

  ! Get the varids of the pressure and temperature netCDF variables.
!  call check( nf90_inq_varid(ncid, PRES_NAME, pres_varid) )
  call check( nf90_inq_varid(ncid, TEMP_NAME, temp_varid) )
  print*, "============ Get the varids of the pressure and temperature netCDF variables."

  ! Read 1 record of NLVLS*NLATS*NLONS values, starting at the beginning 
  ! of the record (the (1, 1, 1, rec) element in the netCDF file).
  count = (/ NLONS, NLATS, NLVLS, 1 /)
  start = (/ 1, 1, 1, 1 /)
  print*, "===Read 1 record of NLVLS*NLATS*NLONS values, starting at the beginning of the record ==="

  ! Read the surface pressure and temperature data from the file, one
  ! record at a time.
  do rec = 1, NRECS
     start(4) = rec
!     call check( nf90_get_var(ncid, pres_varid, pres_in, start = start, &
!                 count = count) )
     call check( nf90_get_var(ncid, temp_varid, temp_in, start, count) )
     
!     i = 0
!     do lvl = 1, NLVLS
!        do lat = 1, NLATS
!           do lon = 1, NLONS
!              if (pres_in(lon, lat, lvl) /= SAMPLE_PRESSURE + i) stop 2
!              if (temp_in(lon, lat, lvl) /= SAMPLE_TEMP + i) stop 2
!              i = i + 1
!           end do
!        end do
!     end do
     ! next record
  end do

      print*, temp_in(30,30,1)

  print*, " Read the surface pressure and temperature data from the file, one record at time."
         
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file.
  call check( nf90_close(ncid) )

  ! If we got this far, everything worked as expected. Yipee! 
  print *,"*** SUCCESS reading example file ", FILE_NAME, "!"

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
end program pres_temp_4D_rd

