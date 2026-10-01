module neutron_math
   
    use precision_variables
    use c_math_bindings ! cbrt call from iso_c_binding
    use philox,    only: philox_fill_uniform_oo ! fills array with (0,1)
    use type_bank, only: material_info, neutron_state, tracking_state, &
                         neutron_tally, triso_dimensions, macro_data

    implicit none

    private
    public :: populate_alpha, gather_pathlengths, spawn_neutrons, &
              event_sampling

    ! module paraeters
    ! philox cbrng key
    integer(id), parameter :: key(2)           = [123_id, 456_id]
    integer(ip), parameter :: material_count   = 3_ip
    ! a = [UO2, C12, SiC] 
    integer(ip), parameter :: a                = [267_ip, 12_ip, 41_ip]
    real(dp),    parameter :: pi               = 3.14159265358979323846_dp
    
    integer(ip), parameter :: u235             = 1_ip
    integer(ip), parameter :: u238             = 2_ip
    
    ! t parameters for watt specturm fission energy sampling
    real(sp),    parameter :: t_u235_thermal   = 1.29_sp
    real(sp),    parameter :: t_u235_fast      = 1.33_sp
    real(sp),    parameter :: t_u238_fast      = 1.40_sp
    real(sp),    parameter :: energy_threshold = 1.0_sp

    ! for event sampling
    integer(ip), parameter :: elastic_scatter = 1_ip
    integer(ip), parameter :: inelastic_scatter = 2_ip
    integer(ip), parameter :: capture = 3_ip
    integer(ip), parameter :: fission = 4_ip
 
    contains

    ! conductor routines that will use elemental sample calls in submodule
    
    subroutine populate_alpha(material, alpha)
        ! in
        type(material_info), intent(in)  :: material
        ! out
        real(sp),            intent(out) :: alpha(3)
        ! local
        real(sp) :: adjust ! adjust from U235 to U238

        adjust = 3.0_sp * (1.0_sp - material%enrichment_fraction)
        
        ! alpha = [UO2(kernel), C12(buffer/pyc), SiC]

        alpha(1) = ((a(1) + adjust - 1.0_sp) / (a(1) + adjust + 1.0_sp)) ** 2
        
        alpha(2) = ((a(2) - 1.0_sp) / (a(2) + 1.0_sp)) ** 2

        alpha(3) = ((a(3) - 1.0_sp) / (a(3) + 1.0_sp)) ** 2
    end subroutine populate_alpha

    pure subroutine scatter_neutrons(neutrons, tracking, randoms)
        ! in
        real(dp),             intent(in)     :: randoms(:)
        ! in out
        type(neutron_state),  intent(in out) :: neutrons
        type(tracking_state), intent(in out) :: tracking
        ! local
        integer(ip) :: i
        real(dp)    :: mu_cm, mu_lab
        ! what will I do with my alpha? 
        
        ! omp flag will go here
        do i = 1_ip, size(neutrons%energy)
            tracking%
        end do
        
    end subroutine scatter_neutrons

    pure subroutine gather_pathlengths(neutrons, tracking, randoms)
        ! in
        real(dp),             intent(in)     :: randoms(:)
        ! in out
        type(neutron_state),  intent(in out) :: neutrons
        type(tracking_state), intent(in out) :: tracking
        ! local
        integer(id)           :: counter(4)
    
        ! the memory for randoms here will be allocated with the size of the
        ! group of neutrons, perhaps when the neutrons split into different
        ! events I might split randoms, or maybe just slice it.

        call philox_fill_uniform_oo(counter, key,    1_id,  &
                                    0.0_dp,  1.0_dp, randoms)

        tracking%current_pathlength =                              &
        distance_to_collision(tracking%current_macro_total, randoms)
    end subroutine gather_pathlengths

    pure subroutine spawn_neutrons(neutrons, material, triso)
        ! in
        type(triso_dimensions), intent(in)  :: triso
        ! in out
        type(neutron_state), intent(in out) :: neutrons
        type(material_info), intent(in out) :: material
        ! local
        integer(id) :: counter(4)
        integer(ip) :: isotope
        integer(ip) :: i
        real(dp)    :: randoms(9)

        counter    = 0_id
        counter(1) = neutrons%generation
    
        do i = 1_ip, size(neutrons%energy)
            counter(2) = i
            call philox_fill_uniform_oo(counter, key,    1_id,  &
                                        0.0_dp,  1.0_dp, randoms)

            call initial_position(randoms(1),          &
                                  randoms(2),          &
                                  randoms(3),          &
                                  triso%kernel_radius, &
                                  neutrons%x(i),       &
                                  neutrons%y(i),       &
                                  neutrons%z(i)        )
            
            call initial_direction(randoms(4),    &
                                   randoms(5),    &
                                   neutrons%u(i), &
                                   neutrons%v(i), &
                                   neutrons%w(i)  )
            
            ! sampling what fissionable isotope we hit
            if (material%enrichment_fraction > randoms(6)) then
                    isotope = u235
            else
                    isotope = u238
            end if

            call initial_energy(randoms(7), randoms(8), randoms(9), &
                                isotope,    neutrons%energy(i)      )
        end do
    end subroutine spawn_neutrons
  
end module neutron_math

