module precision_variables
    ! using C vars for when I add in a C bridge for hardware acceleration
    use, intrinsic :: iso_c_binding, only: c_float, c_double, &
                                           c_int32_t, c_int64_t
    implicit none
    public
    integer, parameter :: sp = c_float   ! single precision real
    integer, parameter :: dp = c_double  ! double precision real
    integer, parameter :: ip = c_int32_t ! single precision integer
    integer, parameter :: id = c_int64_t ! double precision integer

end module precision_variables

