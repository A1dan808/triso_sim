module data_bank

    ! this was written prior to me learning how to do batch allocation
    ! in a convenient way like I have for freeing and initializing some
    ! of the other type arrays

    use precision_variables
    use type_bank, only: xs_data, macro_data, nu_data, auxiliary, performance

    implicit none

    interface check_allocation
        module procedure check_allocation_r1i
        module procedure check_allocation_r1r
        module procedure check_allocation_r1l
    end interface check_allocation
    
    contains

    ! simple memory allocation safety checks
    subroutine check_allocation_r1i(array)
        integer(ip), allocatable, intent(in out) :: array(:)
        if (allocated(array)) deallocate(array)
    end subroutine check_allocation_r1i
    
    subroutine check_allocation_r1r(array)
        real(sp), allocatable, intent(in out) :: array(:)
        if (allocated(array)) deallocate(array)
    end subroutine check_allocation_r1r
    
    subroutine check_allocation_r1l(array)
        logical, allocatable, intent(in out) :: array(:)
        if (allocated(array)) deallocate(array)
    end subroutine check_allocation_r1l
    
    subroutine array_initialization(aux, nu, macro, xs, pref)
        ! in out
        type(auxiliary)       , intent(in out)  :: aux
        type(nu_data)         , intent(in out)  :: nu
        type(macro_data)      , intent(in out)  :: macro
        type(xs_data)         , intent(in out)  :: xs
        type(performance)     , intent(in out)  :: pref

        ! local
        integer(ip) :: i, grid_points
        integer(id) :: temp
        
        ! this is getting reworked very soon
        temp = pref%memory_target_mb * ( 1024 ** 2 )
        pref%num_buckets = temp / (pref%bytes_per_bucket * 2)

        ! all memory safety checks
        call check_allocation(aux%master_grid)
        call check_allocation(aux%hash_map_start)
        call check_allocation(aux%hash_map_end)
        ! U238
        call check_allocation(xs%u238_mt1)
        call check_allocation(xs%u238_mt2)
        call check_allocation(xs%u238_mt4)
        call check_allocation(xs%u238_mt18)
        call check_allocation(xs%u238_mt102)
        call check_allocation(xs%u238_mt452)
        ! U235
        call check_allocation(xs%u235_mt1)
        call check_allocation(xs%u235_mt2)
        call check_allocation(xs%u235_mt4)
        call check_allocation(xs%u235_mt18)
        call check_allocation(xs%u235_mt102)
        call check_allocation(xs%u235_mt452)
        ! C12
        call check_allocation(xs%c12_mt1)
        call check_allocation(xs%c12_mt2)
        call check_allocation(xs%c12_mt102)
        ! O16
        call check_allocation(xs%o16_mt1)
        call check_allocation(xs%o16_mt2)
        call check_allocation(xs%o16_mt102)
        ! Si28
        call check_allocation(xs%si28_mt1)
        call check_allocation(xs%si28_mt2)
        call check_allocation(xs%si28_mt102)
        ! Si29
        call check_allocation(xs%si29_mt1)
        call check_allocation(xs%si29_mt2)
        call check_allocation(xs%si29_mt102)
        ! Si30
        call check_allocation(xs%si30_mt1)
        call check_allocation(xs%si30_mt2)
        call check_allocation(xs%si30_mt102)
        ! Macro
        call check_allocation(macro%kernel_total)
        call check_allocation(macro%buffer_total)
        call check_allocation(macro%pyc_total)
        call check_allocation(macro%sic_total)
        call check_allocation(macro%kernel_scatter)
        call check_allocation(macro%buffer_scatter)
        call check_allocation(macro%pyc_scatter)
        call check_allocation(macro%sic_scatter)
        call check_allocation(macro%kernel_fission)
        ! nu bar
        call check_allocation(nu%u235_mt452)
        call check_allocation(nu%u238_mt452)
        
        ! could be getting reworked very soon here
        grid_points = size(aux%input_union_grid)

        allocate(aux%hash_map_start(pref%num_buckets), source=0_ip)
        allocate(aux%hash_map_end(pref%num_buckets), source=0_ip)
        ! U238
        allocate(xs%u238_mt1(grid_points), source=0.0_sp)
        allocate(xs%u238_mt2(grid_points), source=0.0_sp)
        allocate(xs%u238_mt4(grid_points), source=0.0_sp)
        allocate(xs%u238_mt18(grid_points), source=0.0_sp)
        allocate(xs%u238_mt102(grid_points), source=0.0_sp)
        allocate(xs%u238_mt452(grid_points), source=0.0_sp)
        ! U235
        allocate(xs%u235_mt1(grid_points), source=0.0_sp)
        allocate(xs%u235_mt2(grid_points), source=0.0_sp)
        allocate(xs%u235_mt4(grid_points), source=0.0_sp)
        allocate(xs%u235_mt18(grid_points), source=0.0_sp)
        allocate(xs%u235_mt102(grid_points), source=0.0_sp)
        allocate(xs%u235_mt452(grid_points), source=0.0_sp)
        ! C12
        allocate(xs%c12_mt1(grid_points), source=0.0_sp)
        allocate(xs%c12_mt2(grid_points), source=0.0_sp)
        allocate(xs%c12_mt102(grid_points), source=0.0_sp)
        ! O16
        allocate(xs%o16_mt1(grid_points), source=0.0_sp)
        allocate(xs%o16_mt2(grid_points), source=0.0_sp)
        allocate(xs%o16_mt102(grid_points), source=0.0_sp)
        ! Si28
        allocate(xs%si28_mt1(grid_points), source=0.0_sp)
        allocate(xs%si28_mt2(grid_points), source=0.0_sp)
        allocate(xs%si28_mt102(grid_points), source=0.0_sp)
        ! Si29
        allocate(xs%si29_mt1(grid_points), source=0.0_sp)
        allocate(xs%si29_mt2(grid_points), source=0.0_sp)
        allocate(xs%si29_mt102(grid_points), source=0.0_sp)
        ! Si30
        allocate(xs%si30_mt1(grid_points), source=0.0_sp)
        allocate(xs%si30_mt2(grid_points), source=0.0_sp)
        allocate(xs%si30_mt102(grid_points), source=0.0_sp)
        ! Macro
        ! total
        allocate(macro%kernel_total(grid_points), source=0.0_sp)
        allocate(macro%buffer_total(grid_points), source=0.0_sp)
        allocate(macro%pyc_total(grid_points), source=0.0_sp)
        allocate(macro%sic_total(grid_points), source=0.0_sp)
        ! scatter
        allocate(macro%kernel_scatter(grid_points), source=0.0_sp)
        allocate(macro%buffer_scatter(grid_points), source=0.0_sp)
        allocate(macro%pyc_scatter(grid_points), source=0.0_sp)
        allocate(macro%sic_scatter(grid_points), source=0.0_sp)
        ! fission
        allocate(macro%kernel_fission(grid_points), source=0.0_sp)
    end subroutine array_initialization
end module data_bank

