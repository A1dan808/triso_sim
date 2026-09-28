submodule (type_bank) type_procedures
    implicit none

    contains

    module procedure initialize_neutrons
        call self%free ! initializes generation and current count 
        
        ! I will be optimizing the initialization further
        allocate(self%id(n), source=0_ip)
        allocate(self%energy(n), source=0.0_sp)
        allocate(self%active(n), source=.false.)
        allocate(self%pathlength(n), source=0.0_sp)
        allocate(self%x(n), self%y(n), self%z(n), source=0.0_dp)
        allocate(self%u(n), self%v(n), self%w(n), source=0.0_dp)
    end procedure initialize_neutrons

    module procedure free_neutrons
        if (allocated(self%id))         deallocate(self%id)
        if (allocated(self%energy))     deallocate(self%energy)
        if (allocated(self%pathlength)) deallocate(self%pathlength)
        if (allocated(self%weight))     deallocate(self%weight)
        if (allocated(self%active))     deallocate(self%active)
        if (allocated(self%x))          deallocate(self%x)
        if (allocated(self%y))          deallocate(self%y)
        if (allocated(self%z))          deallocate(self%z)
        if (allocated(self%u))          deallocate(self%u)
        if (allocated(self%v))          deallocate(self%v)
        if (allocated(self%w))          deallocate(self%w)

        self%generation    = 1_ip
        self%current_count = 0_ip
    end procedure free_neutrons

    module procedure initialize_tracking
        call self%free

        allocate(self%region_id(n), source=0_ip)
        allocate(self%material_id(n), source=0_ip)
        allocate(self%current_macro_total(n), source=0.0_sp)
        allocate(self%current_pathlength(n), source=0.0_sp)
    end procedure initialize_tracking

    module procedure free_tracking
        if (allocated(self%region_id)) deallocate(self%region_id)
        if (allocated(self%material_id)) deallocate(self%material_id)
        if (allocated(self%current_macro_total)) &
            deallocate(self%current_macro_total)
        if (allocated(self%current_pathlength)) &
            deallocate(self%current_pathlength)
    end procedure free_tracking
end submodule type_procedures



