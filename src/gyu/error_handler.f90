module error_handler_mod
contains
    subroutine error_mesg(message)
        character(len=*), intent(in) :: message
        print*, "ERROR: ", message
        stop
    end subroutine error_mesg

    subroutine file_check(ios, iomsg, op)
        implicit none
        integer,          intent(in) :: ios
        character(len=*), intent(in) :: iomsg, op
        if (ios == 0) return   ! There was no error, continue
        print*, "Error encountered during " // trim(op)
        print*, "Error code: ", ios
        print*, "Error message: " // trim(iomsg)
        STOP 1
    end subroutine file_check
    
end module error_handler_mod

