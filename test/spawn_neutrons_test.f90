program spawn_neutrons_test
    
    use percision_variables
    use omp_lib
    use neutron_math, only: spawn_neutrons
    use type_bank,    only: neutron_state, material_info, triso_dimentions
    
    implicit none
    
    ! declarations
    integer(ip)            :: n
    type(neutron_state)    :: neutrons
    type(material_info)    :: material
    type(triso_dimentions) :: triso
    
    n = 100000_ip
    triso%kernel_radius = 1.0_sp
    ! initialization
    call neutrons%initialize(n)
    
    call spawn_neutrons(neutrons, material, triso)
    
    print *, 'size sanity check on energy :', size(neutrons%energy)
    print *, 'average initial energy :', real(sum(neutrons%energy) / n, sp)
    print *, 'average x :', real(sum(neutrons%x) / n, sp)
    print *, 'average y :', real(sum(neutrons%y) / n, sp)
    print *, 'average z :', real(sum(neutrons%z) / n, sp)
    print *, 'average u :', real(sum(neutrons%u) / n, sp)
    print *, 'average v :', real(sum(neutrons%v) / n, sp)
    print *, 'average w :', real(sum(neutrons%w) / n, sp)
    print *, 'sum of x y z :', sum(neutrons%x + neutrons%y + neutrons%z)
    print *, 'sum of u v w :', sum(neutrons%u + neutrons%v + neutrons%w)
    ! okay should the sum of each direction component be 0?
    ! why the fuck does the number walk in the negative direction???? 
    ! when sampling x and v the average is lower than the y and z sample
    ! is this a notable statistical issue? would f64 fix this? 
    print *, 'unit vector check on uvw :',     &
    maxval(sqrt((neutrons%u ** 2_ip) +         &
                (neutrons%v ** 2_ip) +         &
                (neutrons%w ** 2_ip))          )
    ! LETS FUCKING GO!!!! unit vetor test passed if valid operations
    ! add more later
    call neutrons%free
end program spawn_neutrons_test





