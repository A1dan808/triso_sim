module type_bank
    use precision_variables
    implicit none
    public
    
    type :: neutron_state
        ! global position and direction
        real(dp),    allocatable :: x(:), y(:), z(:)
        real(dp),    allocatable :: u(:), v(:), w(:)
        real(sp),    allocatable :: energy(:)
        real(sp),    allocatable :: pathlength(:)
        ! auxiliary
        integer(ip)              :: current_count ! write pointer (fission)
        integer(ip)              :: generation
        integer(ip), allocatable :: id(:)
        real(sp),    allocatable :: weight(:)
        logical,     allocatable :: active(:)
        contains
            procedure :: initialize => initialize_neutrons
            procedure :: free       => free_neutrons
    end type neutron_state

    type :: tracking_state
        ! geometry identifiers
        integer(ip), allocatable :: region_id(:)
        integer(ip), allocatable :: material_id(:)
        ! event identifier
        integer(ip), allocatable :: event(:)
        ! sampled value arrays
        real(sp),    allocatable :: current_macro_total(:)
        real(sp),    allocatable :: current_pathlength(:)
        real(dp),    allocatable :: randoms(:) ! for hot loops
        contains
            procedure :: initialize => initialize_tracking
            procedure :: free       => free_tracking
    end type tracking_state

    type :: neutron_tally
        ! will need to add the top three to the namelist config
        integer(ip) :: num_neutrons
        integer(ip) :: num_batches
        integer(ip) :: num_generations
        integer(ip) :: captured
        integer(ip) :: fissioned
        integer(ip) :: leaked 
    end type neutron_tally

    type :: triso_dimensions ! cm
        real(sp) :: kernel_radius
        real(sp) :: buffer_radius
        real(sp) :: inner_pyc_radius
        real(sp) :: sic_radius
        real(sp) :: outer_pyc_radius
    end type triso_dimensions

    type :: geometry_cell
        integer(ip) :: material_id(5)
        integer(ip) :: boundary_count
        integer(ip) :: surface_ids(2)
        ! storing as an f32 should make same kind operations faster
        real(sp)    :: surface_senses(2) ! -1 = inside, 1 = outside
    end type geometry_cell

    type :: geometry_surface
        ! length is the number of layers
        ! In the distant future I will add dynamically configurable layers
        real(sp) :: x_center(5)
        real(sp) :: y_center(5)
        real(sp) :: z_center(5)
        real(sp) :: radius_squared(5)
    end type geometry_surface

    type :: material_info
        ! at% and rho will be inputted via namelist in the .txt file
        real(dp) :: enrichment_fraction
        real(dp) :: rho_uo2
        real(dp) :: rho_buffer
        real(dp) :: rho_pyc
        real(dp) :: rho_sic
    end type material_info

    type :: performance 
        ! performance parameters for the accelerated indexing
        ! resolution of the buckets will be configured base on L2 cache size 
        integer(ip) :: memory_target_mb
        integer(ip) :: bytes_per_bucket
        integer(ip) :: num_buckets
        real(sp)    :: avg_items_per_bucket
        real(sp)    :: max_items_per_bucket
    end type performance

    type :: xs_data
        ! like the above types neutron and tracking state, I will add a free
        ! and an initialize call to the xs data as well
        ! mt1 = total xs
        ! mt2 = elastic scattering xs
        ! mt4 = total inelastic scattering
        ! mt18 = fission
        ! mt102 = radiative capture (abs)
        ! mt452 = nu bar
        ! uranium
        ! U238
        real(sp), allocatable :: u238_mt1(:), u238_mt2(:), u238_mt4(:)
        real(sp), allocatable :: u238_mt18(:), u238_mt102(:), u238_mt452(:)
        ! U235
        real(sp), allocatable :: u235_mt1(:), u235_mt2(:), u235_mt4(:)
        real(sp), allocatable :: u235_mt18(:), u235_mt102(:), u235_mt452(:)
        ! carbon 12
        real(sp), allocatable :: c12_mt1(:), c12_mt2(:), c12_mt102(:)
        ! oxygen 16
        real(sp), allocatable :: o16_mt1(:), o16_mt2(:), o16_mt102(:)
        ! silica
        ! si28
        real(sp), allocatable :: si28_mt1(:), si28_mt2(:), si28_mt102(:)
        ! si29
        real(sp), allocatable :: si29_mt1(:), si29_mt2(:), si29_mt102(:)
        ! si30
        real(sp), allocatable :: si30_mt1(:), si30_mt2(:), si30_mt102(:)
    end type xs_data

    type :: macro_data
        ! Just like the xs type I will add free and initialize routines
        ! and interfaces
        real(sp), allocatable :: kernel_total(:),  kernel_scatter(:)
        real(sp), allocatable :: kernel_fission(:)
        real(sp), allocatable :: buffer_total(:),  buffer_scatter(:)
        real(sp), allocatable :: pyc_total(:),     pyc_scatter(:)
        real(sp), allocatable :: sic_total(:),     sic_scatter(:)
    end type macro_data

    type :: nu_data
        ! considering joining this with the above type
        real(sp), allocatable :: u235_mt452(:)
        real(sp), allocatable :: u238_mt452(:)
    end type nu_data

    type :: auxiliary
        ! for the unionized grids and bucket index accelerator, subject to
        ! change very soon
        real(sp), allocatable :: input_union_grid(:)
        real(sp), allocatable :: master_grid(:)
        integer(ip), allocatable :: hash_map_start(:)
        integer(ip), allocatable :: hash_map_end(:)
    end type
    
    interface
        module subroutine initialize_neutrons(self, n)
            ! in
            integer(ip),          intent(in)     :: n
            ! in out
            class(neutron_state), intent(in out) :: self
        end subroutine initialize_neutrons

        module subroutine free_neutrons(self)
            ! in out
            class(neutron_state), intent(in out) :: self
        end subroutine free_neutrons

        module subroutine initialize_tracking(self, n)
            ! in
            integer(ip),           intent(in)     :: n
            ! in out
            class(tracking_state), intent(in out) :: self
        end subroutine initialize_tracking

        module subroutine free_tracking(self)
            ! in out
            class(tracking_state), intent(in out) :: self
        end subroutine free_tracking
    end interface
end module type_bank

