module data_overhead ! subject to a rename or potential rework
    use precision_variables
    use type_bank

    private

    contains
    ! could be some variable name mismatch errors

    recursive subroutine quicksort(a, first, last)
        ! in
        integer(ip), intent(in)     :: first, last
        ! in out
        real(sp),    intent(in out) :: a(:) ! sp array that needs sorting
        ! local
        real(sp)    :: pivot, temp
        integer(ip) :: i, j

        if (first >= last) return

        ! chose the pivot
        pivot = a((first + last) / 2_ip)
        i = first
        j = last

        ! partition arrays
        do
            do while (a(i) < pivot)
                i = i + 1_ip
            end do

            do while (a(j) > pivot)
                j = j - 1_ip
            end do

            if (i <= j) then
                ! swap arrays
                temp = a(i)
                a(i) = a(j)
                a(j) = temp
                i = i + 1_ip
                j = j - 1_ip
            end if

            if (i > j) exit
        end do

        ! recursively sort the sub arrays
        if (first < j) call quicksort(a, first, j)
        if (i < last)  call quicksort(a, i, last)
    end subroutine quicksort

    subroutine unionized_e_grid(file_list, output_union_grid)
        ! in
        character(len=*), allocatable, intent(in)  :: file_list(:) ! needs work
        real(dp),         allocatable, intent(out) :: output_union_grid(:)
        ! local
        integer(ip) :: grid_points
        integer(ip) :: total_files, i, io_stat, current_count 
        integer(ip) :: total_scratch_size, j, unique_count
        logical               :: valid_file
        ! columns of nuclear data from .dat, extracted
        real(dp)              :: energy, xs
        real(sp), allocatable :: scratch_energies(:)

        total_files = size(file_list)
        total_scratch_size = 500000_ip
        current_count = 0_ip

        allocate(scratch_energies(total_scratch_size), source=0.0_sp)

        do i = 1, total_files
            inquire(file=trim(file_list(i)), exist=valid_file)
            if (.not. valid_file) then
                print *, 'Path to file: "', trim(file_list(i)), &
                         '" is not valid'
                error stop
            end if

            open(unit=10, file=trim(file_list(i)), status='old', &
                 action='read', iostat=io_stat)
            if (io_stat /= 0) then
                print *, 'Error opening file, iostat =', io_stat
                error stop
            end if

            read(10, *, iostat=io_stat)

            ! read file contents and throw energies on scratch pad
            do
                read(10, *, iostat=io_stat) energy, xs
                if (io_stat /= 0) exit ! EOF reached

                current_count = current_count + 1
                if (current_count > total_scratch_size) then
                    print *, 'Scratchpad overflow!', &
                             'Increase total_scratch_size.'
                    error stop
                end if

                scratch_energies(current_count) = energy
            end do
            close(10)
        end do

        ! sort energies after they have been dumped into scratpad
        call quicksort(scratch_energies, 1, current_count)

        ! removing duplicates to create the union set
        unique_count = 0
        if (current_count > 0) then
            unique_count = 1
            do i = 2, current_count
                if (scratch_energies(i) > scratch_energies(unique_count)) then
                    unique_count = unique_count + 1
                    scratch_energies(unique_count) = scratch_energies(i)
                end if
            end do
        end if

        ! allocate output grid to correct size for passing on
        ! this is also subject to being reworked
        allocate(output_union_grid(unique_count))
        output_union_grid = scratch_energies(1:unique_count)

        deallocate(scratch_energies)
    end subroutine unionized_e_grid

    subroutine file_input_map(filename, xs_array, master_grid)
        ! in
        character(len=*), intent(in)     :: filename
        ! in out
        real(sp),         intent(in out) :: master_grid(:)
        real(sp),         intent(in out) :: xs_array(:)
        ! local
        integer(ip) :: grid_points, csv_count, i, k, io_stat
        real(sp)    :: slope
        logical     :: valid_file
        real(sp)    :: energy, xs

        ! dynamic buffers for raw file contents
        real(sp), allocatable :: file_e(:), file_xs(:)

        grid_points = size(master_grid)

        inquire(file=trim(filename), exist=valid_file)
        if (.not. valid_file) then
            print *, 'Error: File not found -> ', trim(filename)
            error stop
        end if

        open(unit=10, file=trim(filename), status='old', action='read', &
             iostat=io_stat)
        if (io_stat /= 0) then
            print *, 'Error opening file: ', trim(filename)
            error stop
        end if

        read(20, *, iostat=io_stat)
        csv_count = 0
        do
            read(20, *, iostat=io_stat) energy, xs
            if (io_stat /= 0) exit ! EOF

            file_count     = file_count + 1_ip
            e(file_count)  = energy
            xs(file_count) = xs
        end do
        close(20)

        if (csv_count == 0_ip) then
            print *, 'Warning: File was empty -> ', trim(filename)
            deallocate(file_e, file_xs)
            return
        end if

        ! two pointer marching interpolation
        k = 1
        do i = 1, grid_points
            ! outside energy range, cross section zero
            if ((master_grid(i) < file_e(1)) .or. &
                (master_grid(i) > file_e(csv_count))) then
                xs_array(i) = 0.0_dp
            else

                do while ((k < file_count - 1_ip) .and. &
                         (file_e(k + 1_ip) <= master_grid(i)))
                    k = k + 1_ip
                end do

                ! check for exact match or calculate slope
                if (master_grid(i) == file_e(k)) then
                    xs_array(i) = file_xs(k)
                else
                    ! linear interpolation
                    slope = (file_xs(k + 1_ip) - file_xs(k)) / &
                            (file_e(k + 1_ip) - file_e(k))

                    xs_array(i) = file_xs(k) + slope         &
                                * (master_grid(i) - file_e(k))
                end if
            end if
        end do

        deallocate(file_e, file_xs)
    end subroutine file_input_map
end module data_overhead

