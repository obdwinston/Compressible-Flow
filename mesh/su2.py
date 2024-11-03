import gmsh

gmsh.initialize()
gmsh.open('mesh/mesh.geo')
gmsh.model.mesh.generate(1)
gmsh.model.mesh.generate(2)
gmsh.option.setNumber('Mesh.Format', 42) # .su2 format
gmsh.option.setNumber('Mesh.SaveAll', 1) # save all elements
gmsh.write('mesh/mesh.su2')
gmsh.finalize()
