module triso_config
    use precision_variables
    use type_bank, only: triso_dimensions, material_info, performance

    ! this module does the input file handling, it might get re written
    ! the check bounds interface and routines might be worth keeping

    implicit none

    private
    public :: map_input_txt

    interface check_bounds
        module procedure check_bounds_real, check_bounds_int
    end interface check_bounds

    contains

    ! bounds checks routines
    impure elemental subroutine check_bounds_real(var_name, val, min_val, &
                                                  max_val)
        character(len=*), intent(in) :: var_name
        real(dp), intent(in) :: val, min_val, max_val
        if (val < min_val .or. val > max_val) then
            print *, 'Invalid value for ','"', trim(var_name),'"', &
                     ' out of allowed range'
            print *, 'allowed range is:', min_val, 'to', max_val
            print *, 'Your current value is:', val, 'fix in config file'
            error stop
        end if
    end subroutine check_bounds_real

    impure elemental subroutine check_bounds_int(var_name, val, min_val, &
                                                 max_val)
        character(len=*), intent(in) :: var_name
        integer(ip), intent(in) :: val, min_val, max_val
        if (val < min_val .or. val > max_val) then
            print *, 'Invalid value for ','"', trim(var_name),'"', &
                     ' out of allowed range'
            print *, 'allowed range is:', min_val, 'to', max_val
            print *, 'Your current value is:', val, 'fix in config file'
            error stop
        end if
    end subroutine check_bounds_int
    
    ! namelist variable importing and kind/bounds checking
    subroutine map_input_txt(filename, triso, material, pref)
        ! in out
        character(len=*)      , intent(in out)  :: filename
        type(triso_dimensions), intent(in out)  :: triso
        type(material_info)   , intent(in out)  :: material
        type(performance)     , intent(in out)  :: pref
        ! local
        character(len=256) :: selected_file
        integer(ip) :: num_args, io_check
        logical :: valid_file
        ! local dummy vars for mapping and bounds checks
        integer(ip) :: num_neutrons
        integer(ip) :: num_generations
        real(dp)    :: kernel_radius
        real(dp)    :: buffer_radius
        real(dp)    :: inner_pyc_radius
        real(dp)    :: sic_radius
        real(dp)    :: outer_pyc_radius
        real(dp)    :: enrichment_fraction
        real(dp)    :: rho_uo2
        real(dp)    :: rho_buffer
        real(dp)    :: rho_pyc
        real(dp)    :: rho_sic
        integer(ip) :: memory_target_mb
        integer(ip) :: bytes_per_bucket

        ! namelist for txt inputs
        namelist /triso_config/ &
        num_neutrons,           &
        num_generations,        &
        kernel_radius,          &
        buffer_radius,          &
        inner_pyc_radius,       &
        sic_radius,             &
        outer_pyc_radius,       &
        enrichment_fraction,    &
        rho_uo2,                &
        rho_buffer,             &
        rho_pyc,                &
        rho_sic,                &
        memory_target_mb,       &
        bytes_per_bucket

        ! check if inputted to CL syntax correctly
        num_args = command_argument_count()
        if (num_args >= 1) then
            call get_command_argument(1, filename)
        else
            print *, 'file not given in CL, try in dir. :', &
                     '"fpm run -- <namelist_var_file.txt>"'
            stop
        end if

        filename = adjustl(filename)
        inquire(file=trim(filename), exist=valid_file)
        if (.not. valid_file) then
            print *, 'File path not valid, try again.'
            stop
        else
            selected_file = filename
            print *, 'File valid, proceeding: ', trim(filename)
        end if

        open(unit=10, file=trim(selected_file), status='old', action='read', &
             iostat=io_check)

        if (io_check /= 0) then
            print *, 'Error opening file, iostat =', io_check
            stop
        end if

        read(10, nml=triso_config, iostat=io_check)

        if (io_check /= 0) then
            print *, 'Error reading namelist in file, iostat =', io_check
            stop
        end if

        close(10)

        ! checking namelist kind/bounds
        call check_bounds('num_neutrons', num_neutrons, 10_ip, 10000000_ip)
        call check_bounds('num_generations', num_generations, 1_ip, 1000_ip)
        call check_bounds('kernel_radius', kernel_radius, 0.0_dp, 20.0_dp)
        call check_bounds('buffer_radius', buffer_radius, 0.0_dp, 20.0_dp)
        call check_bounds('inner_pyc_radius', inner_pyc_radius, 0.0_dp, &
                          20.0_dp)
        call check_bounds('sic_radius', sic_radius, 0.0_dp, 20.0_dp)
        call check_bounds('outer_pyc_radius', outer_pyc_radius, 0.0_dp, &
                          20.0_dp)
        call check_bounds('rho_uo2', rho_uo2, 10.7_dp, 11.0_dp)
        call check_bounds('rho_buffer', rho_buffer, 0.8_dp, 1.2_dp)
        call check_bounds('rho_pyc', rho_pyc, 1.8_dp, 2.30_dp)
        call check_bounds('rho_sic', rho_sic, 3.10_dp, 3.30_dp)
        call check_bounds('enrichment_fraction', enrichment_fraction, 0.0_dp, &
                          1.0_dp)
        call check_bounds('memory_target_mb', memory_target_mb, 1_ip, 100_ip)
        call check_bounds('bytes_per_bucket', bytes_per_bucket, 1_ip, 10_ip)

        ! mapping to global types when bounds are checked
        triso%num_neutrons           = num_neutrons
        triso%num_generations        = num_generations
        triso%kernel_radius          = kernel_radius
        triso%buffer_radius          = buffer_radius
        triso%inner_pyc_radius       = inner_pyc_radius
        triso%sic_radius             = sic_radius
        triso%outer_pyc_radius       = outer_pyc_radius
        material%enrichment_fraction = enrichment_fraction
        material%rho_buffer          = rho_buffer
        material%rho_pyc             = rho_pyc
        material%rho_sic             = rho_sic
        material%rho_uo2             = rho_uo2
        pref%bytes_per_bucket        = bytes_per_bucket
        pref%memory_target_mb        = memory_target_mb

    end subroutine map_input_txt
end module triso_config

