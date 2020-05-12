PROGRAM NR1_READ

  USE netcdf

  IMPLICIT NONE

  INTEGER :: i,j,l,t,ncid,irec

  INTEGER, PARAMETER :: NLATS = 181, NLONS = 360, NLEVS = 37, NTS = 12 !prate

  INTEGER :: var_dimid, var_varid, lev_varid, nlev
  REAL :: var(NLONS, NLATS, NLEVS, NTS), tt(NLEVS-10), lev(NLEVS), ll(NLEVS-10)
  REAL*8 :: offset, factor
  INTEGER :: ndims_in, nvars_in, ngatts_in, unlimdimid_in, lev_id

  INTEGER ::  lev_dimid, varid
  INTEGER :: dimids

  ! User define
  CHARACTER (LEN = *), PARAMETER :: FILE_NAME = "./data/ERA-Interim_W_2018.nc"
  CHARACTER (LEN = *), PARAMETER :: OUTFILE_NAME = "./ERA_Interim_W_AUG.nc"
  CHARACTER (LEN = *), PARAMETER :: LEV_NAME = "level"
  CHARACTER (LEN = *), PARAMETER :: var_name = "w"
  CHARACTER (LEN = *), PARAMETER :: ADD_OFFSET = "add_offset"
  CHARACTER (LEN = *), PARAMETER :: SCALE_FACTOR = "scale_factor"
  CHARACTER (LEN = 20) :: lname, units

  ! Open the file. 
  CALL CHECK( NF90_OPEN(FILE_NAME, nf90_nowrite, ncid) )

  ! netCDF variables, dimensions, and global attributes are in the
  ! file; also the dimension id of the unlimited dimension, if there is one.
  CALL CHECK( NF90_INQUIRE(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )

  ! Get the varids of variable
  CALL CHECK( NF90_INQ_VARID(ncid, LEV_NAME, lev_varid) )
  CALL CHECK( NF90_INQ_DIMID(ncid, LEV_NAME, lev_id) )
  CALL CHECK( NF90_INQUIRE_DIMENSION(ncid, lev_id, len=nlev) )
  CALL CHECK( NF90_INQ_VARID(ncid, var_name, var_varid) )
  
  print*, lev_id, nlev

  ! Read data from the file
  CALL CHECK( NF90_GET_VAR(ncid, lev_varid, lev) )
  CALL CHECK( NF90_GET_VAR(ncid, var_varid, var) )

!   print*, lev(:)
  ll(:) = lev(11:)


  ! Attribute read
  CALL CHECK( NF90_GET_ATT(ncid, var_varid, ADD_OFFSET, offset))
  CALL CHECK( NF90_GET_ATT(ncid, var_varid, SCALE_FACTOR, factor))
  CALL CHECK( NF90_GET_ATT(ncid, var_varid, "units", units))
  CALL CHECK( NF90_GET_ATT(ncid, var_varid, "long_name", lname))


  ! Write data
  Do l = 11, nlev
   tt(l-10) = offset + var(114,60,l,8) * factor
  ENDDO



  ! Close the file
  CALL CHECK( NF90_CLOSE(ncid) )

  PRINT *,"*** SUCCESS reading file! ***"

!============================================================= WRITE NC file

  

 ! Create the netCDF file.
!  CALL CHECK(nf90_create(OUTFILE_NAME,nf90_clobber,ncid))

 ! Define the dimensions.
!  CALL CHECK(nf90_def_dim(ncid, "level",nlev-10,lev_dimid))

 ! Define corrdinate variables
!  CALL CHECK(nf90_def_var(ncid, "lev", nf90_real, lev_dimid, lev_id))


  dimids = lev_dimid
 ! Define variable
!  CALL CHECK(nf90_def_var(ncid, var_name,nf90_float,dimids,varid))

 ! add Units
!  CALL CHECK(nf90_put_att(ncid,lev_id,"units","millibars"))
!  CALL CHECK(nf90_put_att(ncid,lev_id,"long_name","pressure_level"))

!  CALL CHECK(nf90_put_att(ncid,varid,"units",units))
!  CALL CHECK(nf90_put_att(ncid,varid,"long_name",lname))

!  CALL CHECK(nf90_enddef(ncid))  !End DEfinitions

 ! Write Data
!   CALL CHECK(nf90_put_var(ncid, lev_id, ll))
!   CALL CHECK(nf90_put_var(ncid,varid,tt))
!   CALL CHECK(nf90_close(ncid))


CONTAINS
SUBROUTINE CHECK(status)
  INTEGER, INTENT (in) :: status
  IF(status /= NF90_NOERR) THEN
    PRINT *, TRIM(NF90_STRERROR(status))
    STOP "Stopped"
  ENDIF
END SUBROUTINE CHECK
END PROGRAM NR1_READ
