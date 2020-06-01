module file_io_mod
use           global_mod
use         namelist_mod
use    error_handler_mod, only: file_check
contains

    subroutine read_input_data(height_in, temp_in, qv_in, w_in)    !{{{
        real, dimension(:), allocatable, intent(out) :: height_in, temp_in, qv_in, w_in

        integer            :: i, header, nlines, io
        character(len=100) :: iomsg
 
        header = header_data
        nlines = 0 

        open (unit=99, file="INPUT/"//file_name, status='old',    &
              access='sequential', form='formatted', action='read',  &
              iostat=io, iomsg=iomsg)
        call file_check(io, iomsg, "READ")

        do
            read(99,*,iostat=io)
            if (io/=0) exit
            nlines = nlines + 1
        end do

        1 rewind(99)
        allocate(height_in(nlines-header))
        allocate(  temp_in(nlines-header))
        allocate(    qv_in(nlines-header))
        allocate(     w_in(nlines-header))

        if (header >= 1) then
            do i = 1, header
                read(99,*)   ! skip header
            end do
        end if

        if (w_from_data) then
            do i = 1, nlines-header, 1
                read(99,*) height_in(i), temp_in(i), qv_in(i), w_in(i)
            end do
        else
            do i = 1, nlines-header, 1
                read(99,*) height_in(i), temp_in(i), qv_in(i)
            end do
            w_in(:) = w_speed
        end if

        ! Convert units
        ! pressure:     [Pa]     -> [hPa]
        ! temperature:  [C]      -> [K]
        ! mixing ratio: [g kg-1] -> [kg kg-1]
        if ((pres_units .eq. "Pa") .and. (vert_var .eq. "p")) height_in = height_in * 0.01   
        if ( temp_units .eq. "C"   ) temp_in   = temp_in + 273.15
        if (   qv_units .eq. "g/kg") qv_in     = qv_in*0.001   
        ! TODO: processing surface data

        close(99, iostat=io, iomsg=iomsg)
        call file_check(io, iomsg, "CLOSE")
        
    end subroutine read_input_data  !}}}

    subroutine write_data(n)    !{{{
        integer, intent(in) :: n

        integer             :: i
        real, dimension(nt) :: time

        if (n == 1) then 
            ! [m] 1D model full coordinate: [nlev]
            call write_var1d( n, z_full, "z_full.bin" )
            ! [m] 1D model half coordinate: [nlev+1]
            call write_var1d( n, z_half, "z_half.bin" )

            ! [m] droplet radius: [nbin]
            call write_var1d( n,          radius,          "radius.bin" )
            ! [m] radius at boundary: [nbin+1]
            call write_var1d( n, radius_boundary, "radius_boundary.bin" )

            ! [kg] droplet mass: [nbin] x [nlev]
            call write_var2d( n,            mass,            "mass.bin" )
            ! [kg] mass at boundary: [nbin+1] x [nlev]
            call write_var2d( n,   mass_boundary,   "mass_boundary.bin" )
        end if

        if (n == nt) then
            ! [s] Integrated time coordinate: [nt]
            time = (/ (i, i=1, n) /) * dt
            call write_var1d( n,  time,  "time.bin" )
        end if

        ! Interpolated variables
        call write_var1d( n,   Prs,   "Prs.bin" )     ! [hPa] pressure: [nlev]
        call write_var1d( n,     T,     "T.bin" )     ! [K] temperature: [nlev]
        call write_var1d( n, THETA, "THETA.bin" )     ! [K] potential temperature: [nlev]
        call write_var1d( n,    qv,    "qv.bin" )     ! [kg kg] mixing ratio: [nlev]
        call write_var1d( n,     W,     "W.bin" )     ! [m s-1] vertical wind: [nlev+1]

        ! Variables for bins: [nbin] x [nlev]
        call write_var2d( n,    Nr,    "Nr.bin" )     ! [# m-3] number of droplet

        ! Mass tendency: [nbin] x [nlev], [nbin+1] x [nlev]
        call write_var2d( n,  dm_dt,  "dm_dt.bin" )   ! [kg s-1] Mass tendency
        call write_var2d( n, dmb_dt, "dmb_dt.bin" )   ! [kg s-1] at boundary
        
    end subroutine write_data   !}}}

    subroutine write_var1d( n, variable, filename )   !{{{
        integer,            intent(in) :: n
        real, dimension(:), intent(in) :: variable
        character(len=*),   intent(in) :: filename

        ! Output data: Binary format
        if (n == 1) then
            open(21, file="OUTPUT/"//trim(filename), action="write", status="replace")
        else
            open(21, file="OUTPUT/"//trim(filename), action="write", position="append")
        end if

        write(21,*) variable

        close(21)
        
    end subroutine write_var1d    !}}}

    subroutine write_var2d( n, variable, filename )   !{{{
        integer,              intent(in) :: n
        real, dimension(:,:), intent(in) :: variable
        character(len=*),     intent(in) :: filename

        ! Output data: Binary format
        if (n == 1) then
            open(21, file="OUTPUT/"//trim(filename), action="write", status="replace")
        else
            open(21, file="OUTPUT/"//trim(filename), action="write", position="append")
        end if

        write(21,*) variable

        close(21)
        
    end subroutine write_var2d    !}}}

end module file_io_mod
