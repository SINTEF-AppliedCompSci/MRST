import sympy
sympy.init_printing(pretty_print=False)

x, y, z = sympy.symbols('x y z')

dervar = [x, y, z]

c1 = sympy.Rational(1, 2)
c2 = sympy.Rational(2, 3)
u1 = ((x - c1)**2)*(y - c1)**2*(z - c1)**2
u2 = ((x - c1)**2)*(y - c1)**2*(z - c1)**2
u3 = -c2*((x - c1)*(y - c1)**2 + (x - c1)**2*(y - c1))*(z - c1)**3

uvec = list([u1, u2, u3])

gradu = []

# gradu[i][j] gets the derivative of the i-th component with respect to the j-th
# variable
for i in range(3):
    g = []
    g.append(sympy.diff(uvec[0], x))
    g.append(sympy.diff(uvec[1], y))
    g.append(sympy.diff(uvec[2], z))
    gradu.append(g)

# stress = 2*mu*(symmetric gradient) + lambda*trace(symmetric gradient)*I
# The example has been designed such that the trace equal zero.
stress = []
for i in range(3):
    scomp = []
    for j in range(3):
        s = 0.5*(gradu[i][j] + gradu[j][i])
        scomp.append(s)
    stress.append(scomp)

# We compute the divergence
div = []
for j in range(3):
    divcomp = 0
    for i in range(3):
        divcomp = divcomp + sympy.diff(stress[i][j], dervar[i])
    divcomp = sympy.simplify(divcomp)
    div.append(divcomp)
