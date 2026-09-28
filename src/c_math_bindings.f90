module c_math_bindings
    use precision_variables
    implicit none

    ! linking the c cbrt function here with an interface for f32 and f64
    ! is called during initial position sampling

    public :: cbrt

    interface cbrt
        pure real(dp) function c_cbrt_dp(x) bind(c, name="cbrt")
            import :: dp
            implicit none
            real(dp), value :: x
        end function c_cbrt_dp

        pure real(sp) function c_cbrt_sp(x) bind(c, name="cbrtf")
            import :: sp
            implicit none
            real(sp), value :: x
        end function c_cbrt_sp
    end interface
end module c_math_bindings
