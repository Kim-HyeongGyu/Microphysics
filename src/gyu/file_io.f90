module file_io_mod
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

end module file_io_mod
