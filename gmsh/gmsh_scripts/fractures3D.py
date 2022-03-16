from pathlib import Path
import numpy as np
import gmsh


def Graham_scan(x):
    if x.shape == (2, 2):
        dx = np.abs(np.diff(x, axis=0)[0])
        col = np.mod(np.argmin(dx), 2)
        idx = np.argsort(x[:, col])
    else:
        c = np.mean(x, axis=0)
        v = np.arctan2(x[:, 1] - c[1], x[:, 0] - c[0])
        idx = np.argsort(v)
    return x[idx, :]


gmsh.initialize()
model_name = Path(__file__).stem
gmsh.model.add(model_name)
factory = gmsh.model.occ


# Fractures
L = 1000
np.random.seed(0)
num_fractures = 20
min_num_corners = 3
max_num_corners = 10
fractures = []
for k in range(num_fractures):
    num_corners = np.random.randint(min_num_corners, max_num_corners + 1)

    # Random polygon in xy plane at z=L/2
    xy = L / 2 * np.random.random((num_corners, 2)) + L / 4
    xy = Graham_scan(xy)
    p = [factory.addPoint(y[0], y[1], L / 2) for y in xy]
    l = [factory.addLine(p[i - 1], p[i]) for i in range(len(p))]
    cl = factory.addCurveLoop(l)
    s = factory.addPlaneSurface([cl])

    # Translate
    xmean = np.mean(xy, axis=0)
    dx = L / 2 * np.random.random(3) - L / 4
    factory.translate([(2, s)], dx[0], dx[1], dx[2] / 2)
    center = factory.getCenterOfMass(2, s)

    # Rotate
    axis = 2 * np.random.random(3) - 1
    axis /= np.linalg.norm(axis)
    angle = 2 * np.pi * np.random.random()
    factory.rotate(
        [(2, s)], center[0], center[1], center[2], axis[0], axis[1], axis[2], angle
    )

    fractures.append((2, s))
factory.synchronize()

# Matrix
xmin = np.full(3, np.inf)
xmax = np.full(3, -np.inf)
for f in fractures:
    bbox = np.array(factory.getBoundingBox(2, f[1]))
    xmin = np.min([bbox[0:3], xmin], axis=0)
    xmax = np.max([bbox[3:], xmax], axis=0)
L = xmax - xmin
pad = [0.25, 0.25, 0.1]
x0 = [xmin[d] - pad[d] * L[d] for d in range(3)]
dx = [xmax[d] + pad[d] * L[d] - x0[d] for d in range(3)]
matrix = factory.addBox(x0[0], x0[1], x0[2], dx[0], dx[1], dx[2])

# Domain
domain, domain_map = factory.fragment([(3, matrix)], fractures)
factory.synchronize()
matrix = domain_map[0]
fractures = []
for frac in domain_map[1:]:
    if isinstance(frac, list):
        for f in frac:
            fractures.append(f[1])
    else:
        fractures.append(frac[1])
print(f"{num_fractures} fractures split into {len(fractures)} pieces")

# Matrix mesh size
h_max = max(dx) / 10

# Fracture mesh sizes (tangential and normal)
h_t = h_max / 2
h_n = h_max / 5
bump = 10
field = gmsh.model.mesh.field
fields = []
for f in fractures:
    area = factory.getMass(2, f)
    fw = field.add("MathEval")
    fields.append(fw)
    df = field.add("Distance")
    field.setNumbers(df, "SurfacesList", [f])

    num_points = int(max(10, round(area / h_t**2)))
    field.setNumber(df, "Sampling", num_points)

    s = f"{h_max} + ({h_n} - {h_max}) * Exp(-F{df}*F{df} / {bump}^2)"
    field.set_string(fw, "F", s)

# Take the min of all fields
min_field = field.add("Min")
field.setNumbers(min_field, "FieldsList", fields)
field.setAsBackgroundMesh(min_field)

dim = 3
gmsh.model.mesh.generate(dim)

# Set matrix label
matrix_tag = 0
matrix_parts = [dt[1] for dt in matrix]
gmsh.model.addPhysicalGroup(dim, matrix_parts, matrix_tag)

# Set fracture labels as 1, ..., num_fractures
for k in range(1, num_fractures + 1):
    fracture_parts = [dt[1] for dt in domain_map[k]]
    gmsh.model.addPhysicalGroup(dim - 1, fracture_parts, k)

# Write
gmsh.write(model_name + ".msh")
gmsh.write(model_name + ".mesh")
gmsh.write(model_name + ".m")

gmsh.finalize()
