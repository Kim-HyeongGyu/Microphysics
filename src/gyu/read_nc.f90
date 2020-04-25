MODULE read_nc_mod
CONTAINS
  SUBROUTINE read_nc_data(FNAME,vname,l_name,nlev,lev_re,var_re)

  USE netcdf

  IMPLICIT NONE

  INTEGER :: ncid, ll,i
  INTEGER ::  var_varid, lev_varid, lev_id
  INTEGER, intent(in) :: nlev
  INTEGER :: ndims_in, nvars_in, ngatts_in, unlimdimid_in
  real, dimension(nlev) :: var, lev
  real, dimension(nlev), intent(out) :: lev_re,var_re
  CHARACTER(LEN=5),INTENT(IN) :: vname
  CHARACTER(LEN=300) :: FNAME
  CHARACTER(LEN=6) ::  l_name

!============================================================= READ NC file

  ! Open the file. 
  CALL CHECK( NF90_OPEN(trim(FNAME), nf90_nowrite, ncid) )


  ! netCDF variables, dimensions, and global attributes are in the
  ! file; also the dimension id of the unlimited dimension, if there is one.
  CALL CHECK( NF90_INQUIRE(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )
!   print*, ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in

  ! Get the varids of variable
  CALL CHECK( NF90_INQ_VARID(ncid, l_name, lev_varid) )
  CALL CHECK( NF90_INQ_VARID(ncid, trim(vname), var_varid) )


  ! Read data from the file
  CALL CHECK( NF90_GET_VAR(ncid, lev_varid, lev) )
  CALL CHECK( NF90_GET_VAR(ncid, var_varid, var) )


   i = 1
   do ll = nlev,1,-1
    lev_re(i) = lev(ll)
    var_re(i) = var(ll)
     i = i + 1
   end do

  ! Close the file
  CALL CHECK( NF90_CLOSE(ncid) )

  PRINT *,"*** SUCCESS reading file! ***"

  END SUBROUTINE read_nc_data

  ! GRIDDIMS - Get dimensions of a netcdf grid file
  SUBROUTINE griddims(infile,ld_name,nz)
  USE netcdf
  IMPLICIT NONE
  INTEGER, INTENT(OUT)  :: nz
  INTEGER :: ncid, lev_id
  CHARACTER(LEN=300), INTENT(IN) :: infile
  CHARACTER(LEN=6), INTENT(IN) ::  ld_name
  

  ! Open netcdf file
  CALL CHECK (nf90_open(trim(infile), nf90_nowrite, ncid))

  !Inquire about the dimensions
  CALL CHECK( NF90_INQ_DIMID(ncid, ld_name, lev_id) )
  CALL CHECK( NF90_INQUIRE_DIMENSION(ncid,lev_id,len=nz) )


  !Close netcdf file
  CALL CHECK(nf90_close(ncid))

  END SUBROUTINE griddims

  
  SUBROUTINE CHECK(status)
  USE netcdf
  IMPLICIT NONE
  INTEGER, INTENT (in) :: status
  IF(status /= NF90_NOERR) THEN
    PRINT*, TRIM(NF90_STRERROR(status))
    STOP "Stopped"
  ENDIF
  END SUBROUTINE CHECK


END MODULE read_nc_mod


