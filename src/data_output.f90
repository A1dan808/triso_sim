module data_output
    use precision_variables
    use :: type_bank,    only : neutron_state, &
                                neutron_tally, & 
                                triso_dimensions, &
                                material_info, &
                                performance
    contains

    ! this is basically the status='replace' before i knew it existed
    ! this logic was imported from a prior sim, it has had some modification
    ! work done to it however it is far from complete

    subroutine remove_file_if_exists(file_name)
        character(len=*), intent(in) :: file_name
        logical :: file_exists

        file_exists = .false.
        inquire(file=trim(file_name), exist=file_exists)
        ! designed for unix systems as evident by the "rm -f"
        if (file_exists) call execute_command_line('rm -f ' // trim(file_name))
    end subroutine remove_file_if_exists
 
    subroutine output_file(neutron, totals, dimensions, mat_info, pref, &
                           k_eff, file_name)
        ! in
        type(neutron_state), intent(in)     :: neutron
        type(neutron_tally), intent(in)     :: totals
        type(triso_dimensions), intent(in)  :: dimensions
        type(material_info), intent(in)     :: mat_info
        type(performance), intent(in)       :: pref
        real(dp), intent(in)                :: k_eff
        character(len=*), intent(in)        :: file_name
        ! local
        character(len=256) :: date
        character(len=256) :: time
        integer(ip)        :: io_check, f_unit

        ! output file writing
        open(newunit=f_unit, file=trim(file_name), status='replace', &
             iostat=io_check)
        if (io_check /= 0) then
            print *, "Error: Could not open output file ", trim(file_name)
            return
        end if

        call date_and_time(date, time)

        ! header
        write(f_unit, '(A)') &
        '====================================================================='
        write(f_unit, '(A)') '                  TRISO SIMULATOR OUTPUT                              '
        write(f_unit, '(A)') &
        '====================================================================='
        write(f_unit, '(A,A,A,A)') &
        ' Date: ', date(1:4)//'-'//date(5:6)//'-'//date(7:8), &
        ' Time: ', time(1:2)//':'//time(3:4)//':'//time(5:6)
        write(f_unit, '(A)') &
        '---------------------------------------------------------------------'

        ! values from the input file
        write(f_unit, '(A)') ' [1] INPUT PARAMETERS ECHO'
        write(f_unit, '(A)') &
        '---------------------------------------------------------------------'
        write(f_unit, '(A)') ' [a] TRISO DIMENSIONS'
        write(f_unit, '(A)')    ' Kernel Radius (cm):     ', &
        dimensions%kernel_radius
        write(f_unit, '(A)')  ' Buffer Radius (cm):     ', &
        dimensions%buffer_radius
        write(f_unit, '(A)')  ' Inner PyC Radius (cm):  ', &
        dimensions%inner_pyc_radius
        write(f_unit, '(A)')   ' SiC Radius (cm):        ', &
        dimensions%sic_radius
        write(f_unit, '(A)') ' Outer PyC Radius (cm):  ', &
        dimensions%outer_pyc_radius
        write(f_unit, '(A)') ' [b] MATERIAL DENSITY INFORMATION'
        write(f_unit, '(A)')    ' Enrichment Fraction (at):', &
        mat_info%enrichment_fraction
        write(f_unit, '(A)')  ' Kernel Density (g/cm^3): ', &
        mat_info%rho_uo2
        write(f_unit, '(A)')  ' Buffer Density (g/cm^3): ', &
        mat_info%rho_buffer
        write(f_unit, '(A)')  ' PyC Density (g/cm^3):    ', &
        mat_info%rho_pyc
        write(f_unit, '(A)') ' SiC Density (g/cm^3):    ', &
        mat_info%rho_sic
        write(f_unit, '(A)') ' [c] PERFORMANCE PARAMETERS' 
        write(f_unit, '(A)') ' '
        write(f_unit, '(A)') &
        '---------------------------------------------------------------------'
        write(f_unit, '(A)') ' [2] SIMULATION RESULTS'
        write(f_unit, '(A)') &
        '---------------------------------------------------------------------'
        write(f_unit, '(A)') ' '
        write(f_unit, '(A)')  ' Final k_eff:            '
        write(f_unit, '(A)')    ' Input Neutron Count:    ', &
        dimensions%num_neutrons
        write(f_unit, '(A)')  ' Input Generation Count: ', &
        dimensions%num_generations
        write(f_unit, '(A)')  ' Total Absorptions:      ', &
        totals%captured
        write(f_unit, '(A)')  ' Total Fissions:         ', &
        totals%fissioned
        write(f_unit, '(A)')  ' Total Scatters:         '

        write(f_unit, '(A)') ' Total Escapes:          ', &
        totals%leaked
        write(f_unit, '(A)') ' Totals Summed:          '
        write(f_unit, '(A)') ' '
        write(f_unit, '(A)') &
        '====================================================================='
        write(f_unit, '(A)') '                           END OF OUTPUT                              '
        write(f_unit, '(A)') &
        '====================================================================='

        close(f_unit)

    end subroutine output_file
end module data_output

