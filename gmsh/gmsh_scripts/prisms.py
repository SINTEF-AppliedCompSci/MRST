from pathlib import Path
import gmsh

gmsh.initialize()
model_name = Path(__file__).stem
gmsh.model.add(model_name)
factory = gmsh.model.geo

L = 1000
h = L
x0 = 0
y0 = 0
x1 = L
y1 = L
p0 = factory.addPoint(x0, y0, 0, h)
p1 = factory.addPoint(x1, y0, 0, h)
p2 = factory.addPoint(x1, y1, 0, h)
p3 = factory.addPoint(x0, y1, 0, h)
l0 = factory.addLine(p0, p1)
l1 = factory.addLine(p1, p2)
l2 = factory.addLine(p2, p3)
l3 = factory.addLine(p3, p0)
pp = [p0, p1, p2, p3]
ll = [l0, l1, l2, l3]
cl = factory.addCurveLoop(ll)
domain = factory.addPlaneSurface([cl])
factory.extrude([(2, domain)], 0, 0, 100, numElements=[1], recombine=True)
factory.synchronize()

gdim = 3
gmsh.model.mesh.generate(3)
for dim in range(0, gdim + 1):
    for dt in gmsh.model.getEntities(dim):
        gmsh.model.addPhysicalGroup(dim, [dt[1]], dt[1])

gmsh.write(model_name + ".msh")
gmsh.write(model_name + ".m")
gmsh.finalize()
