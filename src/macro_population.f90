module macro_population
    use precision_variables
    use type_bank, only: macro_data, xs_data, material_info

    private
    public :: populate_all_macros

    ! nuclear constants
    ! Avogadro's number, shifted to work with atoms/(barns*cm)
    real(dp), parameter :: na = 0.602214076_dp
    ! atomic masses and wt% from CON 17th edition 
    ! Knolls Atomic Power Laboratory
    real(dp), parameter :: am_u235  = 235.043930_dp
    real(dp), parameter :: am_u238  = 238.050788_dp
    real(dp), parameter :: am_c12   = 12.0000000000_dp
    real(dp), parameter :: am_o16   = 15.9949146196_dp
    real(dp), parameter :: am_si28  = 27.976926532_dp
    real(dp), parameter :: am_si29  = 28.97649470_dp
    real(dp), parameter :: am_si30  = 29.97377017_dp
    ! wt% ratios for silica composition
    real(dp), parameter :: wt_percent_si28 = 92.223_dp
    real(dp), parameter :: wt_percent_si29 = 4.685_dp
    real(dp), parameter :: wt_percent_si30 = 3.092_dp

    contains
    
    subroutine populate_all_macros(macro, xs, material) ! conductor
        ! in
        type(xs_data),       intent(in)  :: xs
        type(material_info), intent(in)  :: material
        ! in out
        type(macro_data), intent(in out) :: macro
        
        call macro_kernel_population(macro, xs, material)
        call macro_buffer_population(macro, xs, material)
        call macro_pyc_population(macro, xs, material)
        call macro_sic_population(macro, xs, material)
    end subroutine populate_all_macros

    subroutine macro_kernel_population(macro, xs, material)
        ! in
        type(xs_data),       intent(in) :: xs
        type(material_info), intent(in) :: material
        ! in out
        type(macro_data), intent(in out) :: macro

        ! local
        real(dp)    :: n_u235, n_u238, n_o16 ! number densities
        real(dp)    :: w                     ! atomic fraction

        w = material%enrichment_fraction
        ! number density calculations
        n_u235 = (material%rho_uo2 * na * w) / &
                 (am_u235 * w + (1.0_dp - w) * am_u238 + 2 * am_o16)

        n_u238 = (material%rho_uo2 * na * (1.0_dp - w)) / &
                 (am_u235 * w + (1.0_dp - w) * am_u238 + 2 * am_o16)

        n_o16  = (material%rho_uo2 * na * 2) / &
                 (am_u235 * w + (1.0_dp - w) * am_u238 + 2 * am_o16)

        macro%kernel_total = (xs%u235_mt1 * n_u235) &
                           + (xs%u238_mt1 * n_u238) &
                           + (xs%o16_mt1  * n_o16)

        macro%kernel_scatter = n_u235 * (xs%u235_mt2 + xs%u235_mt4) &
                             + n_u238 * (xs%u238_mt2 + xs%u238_mt4) &
                             + n_o16  * (xs%o16_mt2) ! will get mt4 later

        macro%kernel_fission = (xs%u235_mt18 * n_u235) &
                             + (xs%u238_mt18 * n_u238)
    end subroutine macro_kernel_population
 
    ! carbon 12 buffer but different density than pyc layer(s)
    subroutine macro_buffer_population(macro, xs, material)
        ! in
        type(xs_data),       intent(in)     :: xs
        type(material_info), intent(in)     :: material
        ! in out
        type(macro_data),    intent(in out) :: macro
        ! local
        real(dp) :: n_buffer

        ! number density calculation 
        n_buffer = (material%rho_buffer * na) / am_c12
        
        macro%buffer_total = n_buffer * xs%c12_mt1

        macro%buffer_scatter = n_buffer * xs%c12_mt2 ! will get mt4 later
    end subroutine macro_buffer_population

    subroutine macro_pyc_population(macro, xs, material)
        ! in
        type(xs_data),       intent(in)     :: xs
        type(material_info), intent(in)     :: material
        ! in out
        type(macro_data),    intent(in out) :: macro
        ! local
        real(dp) :: n_pyc

        ! number density calculation
        n_pyc = (material%rho_pyc * na) / am_c12
        
        macro%pyc_total = n_pyc * xs%c12_mt1

        macro%pyc_scatter = n_pyc * xs%c12_mt2
    end subroutine macro_pyc_population

    subroutine macro_sic_population(macro, xs, material)
        ! in
        type(xs_data),    intent(in)     :: xs
        type(material_info), intent(in)  :: material
        ! in out
        type(macro_data), intent(in out) :: macro
        ! local
        real(dp) :: n_c12, n_si28, n_si29, n_si30

        ! number density calculations
        n_c12  = (material%rho_sic * na) / &
        (am_c12 + am_si28 * wt_percent_si28 + am_si29 * wt_percent_si29 &
        + am_si30 * wt_percent_si30)

        n_si28 = (material%rho_sic * na * wt_percent_si28) / &
        (am_c12 + am_si28 * wt_percent_si28 + am_si29 * wt_percent_si29 &
        + am_si30 * wt_percent_si30)

        n_si29 = (material%rho_sic * na * wt_percent_si29) / &
        (am_c12 + am_si28 * wt_percent_si28 + am_si29 * wt_percent_si29 &
        + am_si30 * wt_percent_si30)

        n_si30 = (material%rho_sic * na * wt_percent_si30) / &
        (am_c12 + am_si28 * wt_percent_si28 + am_si29 * wt_percent_si29 &
        + am_si30 * wt_percent_si30)
        
        macro%sic_total = (n_c12 * xs%c12_mt1)   + (n_si28 + xs%si28_mt1) &
                        + (n_si29 + xs%si29_mt1) + (n_si30 + xs%si30_mt1)

        macro%sic_scatter = (n_c12  * xs%c12_mt2)  + (n_si28 * xs%si28_mt2) &
                          + (n_si29 * xs%si29_mt2) + (n_si30 * xs%si30_mt2)
    end subroutine macro_sic_population
end module macro_population

