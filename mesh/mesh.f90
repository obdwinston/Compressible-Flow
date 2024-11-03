program main
    use iso_fortran_env, only: stdout => output_unit
    use mod_mesh, only: Mesh
    implicit none

    type(Mesh) :: msh
    msh = Mesh()
    
    write(stdout, *) 'Setting mesh...'
    call msh % set_mesh('mesh/mesh.su2', 'mesh/body.txt')
    call msh % save_mesh('mesh/mesh.txt')
    write(stdout, *) 'Done!'
end program main
