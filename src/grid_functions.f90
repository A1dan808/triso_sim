module grid_functions
    use precision_variables
    use type_bank
    
    public

    integer(ip), parameter, private :: mask = int(z'7FFFFFFF', ip)
    
    contains
    
    ! I changed from a log based index call to a bitshift based one
    ! this is one of the most changed modules in the src
    subroutine build_idx_map(master_grid, idx_start, idx_end, &
                              bits_stolen, base_offset        )
        ! in
        real(sp),    intent(in)     :: master_grid(:)
        ! in out
        integer(ip), intent(in out) :: idx_start(:), idx_end(:)
        ! out
        integer(ip), intent(out)    :: bits_stolen, base_offset
        ! local
        integer(ip) :: num_buckets, grid_size
        integer(ip) :: octaves, buckets_per_octaves
        integer(ip) :: i, j, current_bucket, new_bucket, bits_e_min
        integer(id) :: temp
        real(sp)    :: e_min, e_max
        
        ! currently deciding where to put the logic and math for the
        ! idx map arrays
        grid_size = size(master_grid)
        num_buckets = size(idx_start)
        e_min = minval(master_grid)
        e_max = maxval(master_grid)
        
        ! this is grouped with the above statement on needing a rework
        octaves = ceiling(log(e_max / e_min) / log(2.0_sp))
        buckets_per_octaves = num_buckets / octaves
        
        i = 0_ip ! find first non 0 energy value
        do while (master_grid(i) == 0.0_sp)
            ! need to finish the logic here
            i = i + 1_ip
        end do


        if (buckets_per_octaves <= 1_ip) then
            bits_stolen = 0_ip
        else
            ! leadz counts leading zeroes and scales rather than case logic
            bits_stolen = bit_size(buckets_per_octaves) - 1_ip - &
                          leadz(buckets_per_octaves)
        end if

        ! loop to populate the two idx arrays will be written below
    end subroutine build_idx_map

    pure function hash_map_idx(energy, bits_stolen, base_offset) result(idx)
        ! in
        real(sp), intent(in)    :: energy
        integer(ip), intent(in) :: bits_stolen
        integer(ip), intent(in) :: base_offset
        ! return
        integer(ip) :: idx
        ! local
        integer(ip) :: bits

        ! transfer my f32 to an i32, then apply mask
        ! iand is a bitwise and operation
        bits = iand(transfer(energy, bits), mask)

        ! shift it correctly 
        ! odd one because fotran arrays 
        ! subtract base offset
        idx = ishft(bits, bits_stolen - 23_ip) + 1_ip - base_offset
    end function hash_map_idx

    pure function hybrid_search(energy, start_idx, end_idx, master_grid) &
                         result(idx_out)
        ! in
        real(dp),    intent(in) :: energy
        integer(ip), intent(in) :: start_idx, end_idx
        real(dp),    intent(in) :: master_grid(:)
        ! result
        integer(ip)             :: idx_out
        ! local
        integer(ip) :: items, left, right, middle, i

        items = end_idx - start_idx + 1_ip

        ! hybrid search for gpu cores using the L1 cache (not yet tested)
        ! if its under 16 items a linear scan is faster (in theory)
        if (items <= 16) then
        ! linear scan
            i = start_idx
            do while (i < end_idx .and. master_grid(i + 1_ip) <= energy)
                i = i + 1_ip
            end do
            idx_out = i
        else
            ! binary search for larger buckets
            left = start_idx
            right = end_idx
            idx_out = start_idx ! safety fallback, review this later
            do while (left <= right)
                middle = left + (right - left) / 2_ip

                if (master_grid(middle) <= energy) then
                    idx_out = middle
                    left  = middle + 1_ip
                else
                    right = middle - 1_ip
                end if
            end do
        end if
        ! safety catch for max energy in the grid
        if (idx_out == size(master_grid)) idx_out = idx_out - 1_ip
    end function hybrid_search

    pure function lin_interp(energy, lower_idx, chosen_grid, master_grid) &
                      result(y)
        ! dedicated linear interpolation function subject to change, could
        ! compute slope arrays as overhead unless memory forbids it
        ! in
        integer(ip), intent(in) :: lower_idx ! idx_out from above subroutine
        real(sp),    intent(in) :: chosen_grid(:), master_grid(:), energy
        ! result
        real(sp)                :: y
        ! local
        real(sp)            :: x_1, y_1, x_2, y_2, x
        real(dp), parameter :: tol = 1e-14_dp

        x_1 = master_grid(lower_index)
        x_2 = master_grid(lower_index + 1_ip)
        y_1 = chosen_grid(lower_index)
        y_2 = chosen_grid(lower_index + 1_ip)
        x   = energy
        
        ! inspect further
        if ((abs(x_1) < tol) .or. (abs(x_2) < tol)) then
            y = 0.0_sp
        else
            y = y_1 + ((y_2 - y_1)/(x_2 - x_1)) * (x - x_1)
        end if
    end function lin_interp

end module grid_functions
