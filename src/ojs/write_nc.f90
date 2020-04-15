MODULE write_nc_mod
CONTAINS

  SUBROUTINE write_nc_data(ONAME,vname,nt,nz,z_full,var)

  USE netcdf

  IMPLICIT NONE

  INTEGER :: ncid, rec
  INTEGER :: varid, lev_dimid, lev_id, time_dimid, time_id
  INTEGER, intent(in) :: nz, nt
  real, dimension(nz),intent(in) :: z_full
  real, dimension(nz,nt),intent(in) :: var
  CHARACTER(LEN=5),INTENT(IN) :: vname
  CHARACTER(LEN=300),INTENT(IN) :: ONAME

  INTEGER, dimension(2) :: dimids, start, count
  CHARACTER (LEN = 30) :: lname, units


!============================================================= WRITE NC file

  

 ! Create the netCDF file.
  CALL CHECK(nf90_create(trim(ONAME),nf90_clobber,ncid))

 ! Define the dimensions.
  CALL CHECK(nf90_def_dim(ncid, "level",nz,lev_dimid))
  CALL CHECK(nf90_def_dim(ncid, "time",NF90_UNLIMITED,time_dimid))

 ! Define corrdinate variables
  CALL CHECK(nf90_def_var(ncid, "hgt", nf90_real, lev_dimid, lev_id))
  CALL CHECK(nf90_def_var(ncid, "time", nf90_real, time_dimid, time_id))


  dimids = (/lev_dimid, time_dimid/)
 ! Define variable
  CALL CHECK(nf90_def_var(ncid,vname,nf90_float,dimids,varid))

 ! add Units
  CALL CHECK(nf90_put_att(ncid,lev_id,"units","m"))
  CALL CHECK(nf90_put_att(ncid,lev_id,"long_name","height[m]"))

  if (vname == "T") then
      units = "K"
      lname = "Temperature"
  else if (vname == "q") then
      units = "kg kg**-1"
      lname = "Specific humidity"
  else 
      print*, "Error!! Check the variable name !!"
  endif

  CALL CHECK(nf90_put_att(ncid,varid,"units",units))
  CALL CHECK(nf90_put_att(ncid,varid,"long_name",lname))


  CALL CHECK(nf90_enddef(ncid))  !End DEfinitions

!   count = (/nz,1/)
!   start = (/1,1/)
   
 ! Write Data
   CALL CHECK(nf90_put_var(ncid, lev_id, z_full))
!   do rec = 1, nt
!   start(2) = rec
   CALL CHECK(nf90_put_var(ncid,varid,var))
!   CALL CHECK(nf90_put_var(ncid,varid,var,start=start,count= count))
!   end do
   CALL CHECK(nf90_close(ncid))

  PRINT *,"*** SUCCESS writing file! ***"
  END SUBROUTINE write_nc_data

SUBROUTINE CHECK(status)
   USE netcdf
   IMPLICIT NONE
   INTEGER, INTENT (in) :: status
  IF(status /= NF90_NOERR) THEN
    PRINT *, TRIM(NF90_STRERROR(status))
    STOP "Stopped"
  ENDIF
END SUBROUTINE CHECK
END MODULE write_nc_mod
