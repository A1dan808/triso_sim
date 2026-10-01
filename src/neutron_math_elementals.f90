submodule (neutron_math) neutron_math_elementals
    implicit none

contains

! elemental sampling routines and functions

elemental subroutine initial_position(random_1, random_2, random_3, &
                                          radius, x, y, z)
        ! in
        real(sp), intent(in)     :: radius
        real(dp), intent(in)     :: random_1, random_2, random_3
        ! in out
        real(dp), intent(in out) :: x, y, z
        ! local
        real(dp)                 :: r, cos_theta, sin_theta, phi
                
        r = radius * cbrt(random_1)

        cos_theta = 2.0_dp * random_2 - 1.0_dp
        sin_theta = sqrt(1.0_dp - cos_theta ** 2_ip) ! make safer?

        phi = 2.0_sp * pi * random_3

        x = r * sin_theta * cos(phi)
        y = r * sin_theta * sin(phi)
        z = r * cos_theta
    end subroutine initial_position

    elemental subroutine initial_direction(random_1, random_2, u, v, w)
        ! in
        real(dp), intent(in)     :: random_1, random_2
        ! in out
        real(dp), intent(in out) :: u, v, w
        ! local
        real(dp)                 :: phi, sin_theta, cos_theta

        w = 2.0_dp * random_1 - 1.0_dp 

        sin_theta = sqrt(1.0_dp - w**2_ip)

        phi = 2.0_dp * pi * random_2
        
        u = sin_theta * cos(phi)
        v = sin_theta * sin(phi)
    end subroutine initial_direction

    elemental subroutine initial_energy(random_1, random_2, random_3,   &
                                        isotope, energy, incident_energy)
        ! in
        integer(ip), intent(in)        :: isotope
        real(dp),    intent(in)        :: random_1, random_2, random_3
        real(sp), optional, intent(in) :: incident_energy
        ! in out
        ! out
        real(sp), intent(out)          :: energy
        ! local
        real(sp)                       :: t ! T

        ! select case and if statements for sampling what isotope we choose
        ! this will use one of the Ts from the parent module
        ! the select case and if system will be improved upon later when
        ! I'm optimizing

        select case (isotope)
        case(u235)
            if(present(incident_energy)) then
                if (incident_energy > energy_threshold) then
                    t = t_u235_fast
                else
                    t = t_u235_thermal
                end if
            else
                t = t_u235_thermal
            end if
        case(u238)
            ! no thermal u238
            t = t_u238_fast
        case default
            t = t_u235_thermal
        end select
 
        energy = t * (-log(random_1) - log(random_2) &
               * cos(0.5_dp * pi * random_3) ** 2_ip)
    end subroutine initial_energy
   
    elemental function distance_to_collision(macro_xs, random) result(d_c)
        ! in
        real(sp), intent(in) :: macro_xs
        real(dp), intent(in) :: random
        ! result
        real(sp) :: d_c

        d_c = real((-(log(random)) / macro_xs), sp)
    end function distance_to_collision

    elemental function new_weight(old_weight, macro_scatter, macro_total &
                                  macro_fission)            result(weight)
        ! for implicit capture, varience reduction using weights
        ! in
        real(sp),           intent(in) :: old_weight 
        real(sp),           intent(in) :: macro_scatter, macro_total
        real(sp), optional, intent(in) :: macro_fission
        ! result
        real(sp) :: weight

        if (present(macro_fission)) then
            weight = old_weight * ((macro_scatter + macro_fission) &
                                   / macro_total)
        else
            weight = old_weight * (macro_scatter / macro_total)
        end if
    end function new_weight
    
    elemental subroutine elastic_kinematics(energy, random, a, alpha, mu_lab)
        ! incomplete, will be reworked stopping here for time being.
        ! in
        integer(ip), intent(in) :: a
        real(dp),    intent(in) :: random
        real(sp),    intent(in) :: alpha
        ! in out
        real(sp), intent(in out) :: energy
        ! local
        real(sp) :: cos_mew_cm

        cos_mew_cm = real((1.0_sp - 2.0_sp * random), sp)

        new_energy = 0.5_sp * old_energy                             &
                   * ((1.0_sp - alpha) *  cos_mew_cm + 1.0_sp + alpha)
    end subroutine elastic_kinematics
    
    elemental function fission_births(nu_bar, random) result(number_births)
        ! will return an integer amount of neutrons to spawn from fission
        ! in
        real(sp), intent(in) :: nu_bar
        real(dp), intent(in) :: random
        ! result
        integer(ip)          :: number_births
        ! local
        integer(ip)          :: floor_integer, rounding_int
        real(sp)             :: mantissa
        real(sp), parameter  :: tolerance = 1.0e-6_sp

        floor_integer = floor(nu_bar)

        mantissa = nu_bar - floor_integer

        if (mantissa < tolerance) then
            number_births = floor_integer
            return
        end if
        
        if (mantissa > random) then
            number_births = floor_integer + 1_ip
        else
            number_births = floor_integer
        end if
    end function fission_births
    
    elemental function event_cdf(random, macro_total, macro_inelastic, &
                                 macro_elastic, macro_capture,         &
                                 macro_fission)            result(event)
        ! in
        real(dp),           intent(in) :: random
        real(sp),           intent(in) :: macro_total
        real(sp),           intent(in) :: macro_elastic
        real(sp),           intent(in) :: macro_inelastic
        real(sp),           intent(in) :: macro_capture
        real(sp), optional, intent(in) :: macro_fission
        ! result
        integer(ip) :: event
        ! local
        real(sp) :: scaled_random, cumulative_xs

        scaled_random = real(random, sp) * macro_total
        
        ! now with faster flat CDF
        ! check elastic
        cumulative_xs = macro_elastic
        if (scaled_random < cumulative_xs) then
            event = elastic_scatter
            return
        end if

        ! check inelastic
        cumulative_xs = cumulative_xs + macro_inelastic
        if (scaled_random < cumulative_xs) then
            event = inelastic_scatter
            return
        end if

        ! check capture
        cumulative_xs = cumulative_xs + macro_capture
        if (scaled_random < cumulative_xs) then
            event = capture
            return
        end if

        ! check fission if applicable
        if (present(macro_fission)) then
            cumulative_xs = cumulative_xs + macro_fission
            if (scaled_random < cumulative_xs) then
                event = fission
            end if
        end if

    end function event_cdf

    elemental subroutine scatter_direction(u, v, w, mu_lab, random)
        ! in
        real(dp), intent(in)     :: random
        real(dp), intent(in)     :: mu_lab
        ! in out
        real(dp), intent(in out) :: u, v, w
        ! local, temps to reduce cycles
        real(dp) :: phi, mu_sqrt, w_sqrt, cos_phi, sin_phi, norm
        
        ! sample azimuthal angle (0,2pi)
        phi = random * 2.0_dp * pi

        ! calculate and map temps (safe square roots)
        w_sqrt  = sqrt(max(0.0_dp, 1.0_dp - w**2_ip))
        mu_sqrt = sqrt(max(0.0_dp, 1.0_dp - mu_lab**2_ip))
        sin_phi = sin(phi)
        cos_phi = cos(phi)

        ! neutron is not on the z axis
        if (abs(w) < 0.999999_dp) then
            
            u = u*mu_lab + (mu_sqrt * (u*w*cos_phi - v*sin_phi)) &
              / w_sqrt

            v = v*mu_lab + (mu_sqrt * (v*w*cos_phi + u*sin_phi)) &
              / w_sqrt
            
            w = w*mu_lab - mu_sqrt * w_sqrt * cos_phi

        else ! neutron is on the z axis
            
            u = mu_sqrt * cos_phi

            v = mu_sqrt * sin_phi

            w = sign(1.0_dp, w) * mu_lab

        end if

        ! re normalize direction vector to prevent drift, optimize later
        norm = 1.0_dp / sqrt(u**2_ip + v**2_ip + w**2_ip)
        u = u * norm
        v = v * norm
        w = w * norm

    end subroutine scatter_direction

end submodule neutron_math_elementals
