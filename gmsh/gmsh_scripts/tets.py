from pathlib import Path
import gmsh

gmsh.initialize()
model_name = Path(__file__).stem
gmsh.model.add(model_name)
factory = gmsh.model.occ

p0 = factory.addPoint(0, 0, 0)
p1 = factory.addPoint(1, 0, 0)
p2 = factory.addPoint(0, 1, 0)
p3 = factory.addPoint(0, 1, 1)
l0 = factory.addLine(p0, p1)
l1 = factory.addLine(p0, p2)
l2 = factory.addLine(p1, p2)
l3 = factory.addLine(p0, p3)
l4 = factory.addLine(p1, p3)
l5 = factory.addLine(p2, p3)
cl0 = factory.addCurveLoop([l0, l1, l2])
cl1 = factory.addCurveLoop([l0, l3, l4])
cl2 = factory.addCurveLoop([l1, l3, l5])
cl3 = factory.addCurveLoop([l2, l4, l5])
sf = [factory.addSurfaceFilling(cl) for cl in [cl0, cl1, cl2, cl3]]
sl = factory.addSurfaceLoop(sf)
factory.addVolume([sl])

factory.synchronize()

gmsh.option.setNumber("Mesh.MeshSizeFactor", 5)

gdim = 3
gmsh.model.mesh.generate(gdim)

for dim in range(0, gdim + 1):
    for dt in gmsh.model.getEntities(dim):
        gmsh.model.addPhysicalGroup(dim, [dt[1]], dt[1])

gmsh.write(model_name + ".msh")
gmsh.write(model_name + ".m")
gmsh.finalize()
