module file_io_mod
use            global_mod
implicit none
    
contains

    subroutine read_namelist()  ! {{{
        ! Set variables
        namelist /main_nml/     t_final
        namelist /dynamics_nml/ num_levels, top_of_atmosphere,  &
                                vertical_grid, vertical_advect, &
                                CFL_condition
        namelist /physics_nml/  rmin, rratio, nbin, Nc, qc,     &
                                dist_type

        open  (unit = 8, file = 'input.nml', delim = 'apostrophe')
        read  (unit = 8, nml  = main_nml) 
        read  (unit = 8, nml  = dynamics_nml) 
        read  (unit = 8, nml  = physics_nml) 
        close (unit = 8)

        ! Overide default values for optional arguments
        nz      = num_levels
        ztop    = top_of_atmosphere   ! [m]
    end subroutine read_namelist    ! }}}


    subroutine read_data_init() ! {{{
        allocate(Tinit(nz), w(nz), qinit(nz))
        ! Temperature   [K]
        Tinit(:) = (/ (I+273, I = nz,1,-1) /)  ! lapse rate 1K/km

        ! Vertical wind [m s-1]
        w      = sin( (/ (I, I = 1,nz*2,2) /) / 10. )
        w      = w*0. + 2.

        ! Mixing ratio  [kg kg-1]
        qinit(:) = w*0.
        qinit(5) = 10.

        ! initialization
        ! open  (91,file='./INPUT/T.dat',access='direct',recl=...)
        ! read  (91,rec=1) Tinit
        ! close (91)
        ! open  (92,file='./INPUT/q.dat',access='direct',recl=...)
        ! read  (92,rec=1) qinit
        ! close (92)
        ! open  (92,file='./INPUT/w.dat',access='direct',recl=...)
        ! read  (92,rec=1) w
        ! close (92)
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
