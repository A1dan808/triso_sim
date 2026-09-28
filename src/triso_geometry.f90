module triso_geometry
    use precision_variables
    use type_bank, only: geometry_surface, geometry_cell, triso_dimensions

    ! using basic flat CSG, flat for gpu parallelism
    implicit none
    
    private
    public :: build_csg, populate_surface, populate_cells, &
                  evaluate_surface 

    ! hardcoded configuration for now, will be dynamic in the future
    integer(ip), parameter :: kernel    = 1_ip
    integer(ip), parameter :: buffer    = 2_ip
    integer(ip), parameter :: inner_pyc = 3_ip
    integer(ip), parameter :: sic       = 4_ip
    integer(ip), parameter :: outer_pyc = 5_ip

    integer(ip) :: layers = [kernel, buffer, inner_pyc, sic, outer_pyc]

    contains
    
    subroutine build_csg(cell, surface, dimensions)
        ! in
        type(triso_dimensions), intent(in)     :: dimensions 
        ! in out
        type(geometry_surface), intent(in out) :: surface
        type(geometry_cell),    intent(in out) :: cell(size(layers))
        ! local
        integer(ip) :: i
        real(sp)    :: radius

        do i = 1_ip, size(layers)
            select case(layers(i))
                case(kernel)
                    radius = dimensions%kernel_radius
                case(buffer)
                    radius = dimensions%buffer_radius
                case(inner_pyc)
                    radius = dimensions%inner_pyc_radius
                case(sic)
                    radius = dimensions%sic_radius
                case(outer_pyc)
                    radius = dimensions%outer_pyc_radius
            end select
            
            ! populate values for math surface and logical volume rules
            call populate_surface(surface, 0.0_sp, 0.0_sp, 0.0_sp, i, radius)
            call populate_cells(cell(i), i)
        end do
    end subroutine build_csg

    pure subroutine populate_surface(surface, x, y, z, layer, layer_radius)
        ! in
        integer(ip),            intent(in)     :: layer
        real(sp),               intent(in)     :: layer_radius
        real(sp),               intent(in)     :: x, y, z
        ! in out
        type(geometry_surface), intent(in out) :: surface
        
        ! assign global center positions
        surface%x_center(layer) = x
        surface%y_center(layer) = y
        surface%z_center(layer) = z
        surface%radius_squared(layer) = layer_radius ** 2_ip
    end subroutine populate_surface

    pure subroutine populate_cells(cell, layer)
        ! in
        integer(ip),         intent(in)     :: layer
        ! in out
        type(geometry_cell), intent(in out) :: cell
       
        if (layer == 1_ip) then
            ! for the most inner layer (kernel)
            cell%boundary_count    = 1_ip
            cell%surface_ids(1)    = 1_ip
            cell%surface_senses(1) = -1.0_sp
        else 
            ! all other layers
            cell%boundary_count    = 2_ip
            ! bound 1, inner
            cell%surface_ids(1)    = layer - 1_ip
            cell%surface_senses(1) = 1.0_sp
            ! bound 2, outer
            cell%surface_ids(2)    = layer
            cell%surface_senses(2) = -1.0_sp
        end if 
    end subroutine populate_cells

    elemental function evaluate_surface(x, y, z, x_c, y_c, z_c,  &
                                        radius_squared)          &
                                        result(surface_evaluation)
        ! in
        real(dp), intent(in) :: x, y, z
        real(sp), intent(in) :: x_c, y_c, z_c, ! centers
        real(sp), intent(in) :: radius_squared
        ! result
        real(sp) :: surface_evaluation ! [-1, 1]
        
        ! f(x,y,z) = (x-x_c)^2 + (y-y_c)^2 + (z-z_c)^2 - R^2
        surface_evaluation = real((x-x_c)**2_ip + (y-y_c)**2_ip + (z-z_c)**2_ip 
                           - radius_squared, sp) ! is this the best fix? 
    end function evaluate_surface 
    
end module triso_geometry

