module transport
    use precision_variables
    use type_bank, only: neutron_state, triso_dimensions, material_info, &
                         neutron_tally
    use neutron_math, only: spawn_neutron, spawn_fission_neutron
    implicit none

    private
    public :: transport
    
    contains

    subroutine transport(neutron, bank, config, material, tally)
        ! transport is under construction, will be event based
        ! in
        type(triso_dimensions), intent(in)     :: config
        type(material_info),    intent(in)     :: material
        ! in out
        type(neutron_state),    intent(in out) :: neutron, bank
        type(neutron_tally),    intent(in out) :: tally
        ! local
        integer(ip) :: i, j
        
        ! as stated else where batches and generations will be in namelist .txt
        do i = 1_ip, tally%num_batches 
            do j = 1_ip, tally%num_generations
            end do
        end do
    end subroutine transport
end module transport

