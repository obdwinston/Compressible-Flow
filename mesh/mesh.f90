program main
    use mod_mesh, only: Mesh
    implicit none

    type(Mesh) :: msh
    msh = Mesh()
    
    call msh % set_mesh('mesh/mesh.su2', 'mesh/body.txt')
    call msh % save_mesh('mesh/mesh.txt')
end program main
