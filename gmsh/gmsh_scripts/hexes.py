from pathlib import Path
import gmsh

gmsh.initialize()
model_name = Path(__file__).stem
gmsh.model.add(model_name)
factory = gmsh.model.occ

L = 10
box = factory.addBox(0, 0, 0, L, L, L)
factory.synchronize()

num_nodes = 3
for c in gmsh.model.getEntities(1):
    gmsh.model.mesh.setTransfiniteCurve(c[1], num_nodes)
for s in gmsh.model.getEntities(2):
    gmsh.model.mesh.setTransfiniteSurface(s[1])
    gmsh.model.mesh.setRecombine(s[0], s[1])
    gmsh.model.mesh.setSmoothing(s[0], s[1], 100)
gmsh.model.mesh.setTransfiniteVolume(box)
gmsh.model.mesh.setRecombine(3, box)

gdim = 3
gmsh.model.mesh.generate(gdim)

for dim in range(0, gdim + 1):
    for dt in gmsh.model.getEntities(dim):
        gmsh.model.addPhysicalGroup(dim, [dt[1]], dt[1])

gmsh.write(model_name + ".msh")
gmsh.write(model_name + ".m")
gmsh.finalize()
