module model_initialize_mod
use              global_mod
use            namelist_mod
use             file_io_mod, only: read_input_data
use     vert_coordinate_mod, only: compute_vert_coord, &
                                   interpolate_1d
use        microphysics_mod, only: initialize_microphysics
contains

    subroutine model_initialize()

        call read_namelist()

        ! Read data from INPUT/...
        call read_input_data( height_in, temp_in, qv_in, w_in )

        ! Setup vertical info
        call compute_vert_coord( ztop, zbottom, nz, grid_type, &
                                 z_full, z_half, dz)
        call interpolate_1d( height_in, temp_in, qv_in, w_in,  &
                             z_full, z_half, Prs, T, THETA,    &
                             qv, W, THETA_sfc, qv_sfc          )

        ! Setup droplet distribuion
        call initialize_microphysics( radius, radius_boundary, & 
                                      mass, mass_boundary,     & 
                                      Nr, dm_dt, dmb_dt        )   

        call write_model_setup_info()

        call model_initialize_end()

    end subroutine model_initialize

    subroutine model_initialize_end()   !{{{

        ! Deallocate input data
        if (allocated(height_in)) deallocate(height_in)
        if (allocated(  temp_in)) deallocate(  temp_in)
        if (allocated(    qv_in)) deallocate(    qv_in)
        if (allocated(     w_in)) deallocate(     w_in)

    end subroutine model_initialize_end !}}}

    subroutine model_close()    !{{{

        ! Deallocate global variables
        if (allocated(         z_full)) deallocate(         z_full)
        if (allocated(         z_half)) deallocate(         z_half)

        if (allocated(         radius)) deallocate(         radius)
        if (allocated(radius_boundary)) deallocate(radius_boundary)
        if (allocated(           mass)) deallocate(           mass)
        if (allocated(  mass_boundary)) deallocate(  mass_boundary)
        if (allocated(             Nr)) deallocate(             Nr)

        ! Deallocate advecttion variables
        ! if (allocated(           zwts)) deallocate(           zwts)
        ! if (allocated(            dzs)) deallocate(            dzs)

        print*, "Successfully run!"
        
    end subroutine model_close  !}}}

end module model_initialize_mod
