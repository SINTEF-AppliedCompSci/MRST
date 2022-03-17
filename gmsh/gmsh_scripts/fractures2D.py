from pathlib import Path
import numpy as np
import gmsh

gmsh.initialize()
model_name = Path(__file__).stem
gmsh.model.add(model_name)
factory = gmsh.model.occ

# Matrix
L = 1000
matrix = factory.addRectangle(0, 0, 0, L, L)

# Fractures
num_fractures = 20
np.random.seed(0)
x = L * np.random.random((2 * num_fractures, 2))
fractures = []
for k in range(num_fractures):
    p0 = factory.addPoint(x[2 * k, 0], x[2 * k, 1], 0)
    p1 = factory.addPoint(x[2 * k + 1, 0], x[2 * k + 1, 1], 0)
    l = factory.addLine(p0, p1)
    fractures.append((1, l))
factory.synchronize()

# Domain
domain, domain_map = factory.fragment([(2, matrix)], fractures)
factory.synchronize()
fractures = []
for frac in domain_map[1:]:
    if isinstance(frac, list):
        for f in frac:
            fractures.append(f[1])
    else:
        fractures.append(frac[1])
print(f"{num_fractures} fractures split into {len(fractures)} pieces")

# Matrix mesh size
h_max = L / 10

# Fracture mesh sizes (tangential and normal)
h_t = L / 50
h_n = L / 100
bump = 50
field = gmsh.model.mesh.field
fields = []
for f in fractures:
    length = factory.getMass(1, f)
    fw = field.add("MathEval")
    fields.append(fw)
    df = field.add("Distance")
    field.setNumbers(df, "CurvesList", [f])

    num_points = int(max(2, round(length / h_t)))
    field.setNumber(df, "Sampling", num_points)

    s = f"{h_max} + ({h_n} - {h_max}) * Exp(-F{df}*F{df} / {bump}^2)"
    field.set_string(fw, "F", s)

# Take the min of all fields
min_field = field.add("Min")
field.setNumbers(min_field, "FieldsList", fields)
field.setAsBackgroundMesh(min_field)

# Enforce mesh size from field only
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

# Quad mesh
gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 0)
gmsh.option.setNumber("Mesh.RecombineAll", 1)

# Alternative quad mesh by triangle subdivision
# gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)

dim = 2
gmsh.model.mesh.generate(dim)

# Set matrix label
matrix_tag = 0
matrix_parts = [dt[1] for dt in domain_map[0]]
gmsh.model.addPhysicalGroup(dim, matrix_parts, matrix_tag)

# Set fracture labels as 1, ..., num_fractures
for k in range(1, num_fractures + 1):
    fracture_parts = [dt[1] for dt in domain_map[k]]
    gmsh.model.addPhysicalGroup(dim - 1, fracture_parts, k)

# Write
gmsh.write(model_name + ".msh")
gmsh.write(model_name + ".m")

gmsh.finalize()
