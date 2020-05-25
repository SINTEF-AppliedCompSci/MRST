# compute external force for given displacement function

import sympy
sympy.init_printing(pretty_print=False)

x, y, z = sympy.symbols('x y z')
mu, lamb = sympy.symbols('mu0 lamb0')

dervar = [x, y, z]

c1 = sympy.Rational(1, 2)
c2 = sympy.Rational(2, 3)

d = 3
if d == 3:
    u1 = ((x - c1)**2)*(y - c1)**2*(z - c1)**2
    u2 = ((x - c1)**2)*(y - c1)**2*(z - c1)**2
    u3 = -c2*((x - c1)*(y - c1)**2 + (x - c1)**2*(y - c1))*(z - c1)**3
    uvec = list([u1, u2, u3])
elif d == 2:
    u1 = (x - c1)**2*(y - c1)**2
    u2 = -c2*(x - c1)*(y - c1)**3
    uvec = list([u1, u2])


gradu = []

# gradu[i][j] gets the derivative of the i-th component with respect to the j-th
# variable
for i in range(d):
    g = []
    for j in range(d):
        g.append(sympy.diff(uvec[i], dervar[j]))
    gradu.append(g)

# stress = 2*mu*(symmetric gradient) + lambda*trace(symmetric gradient)*I
# The example has been designed such that the trace equal zero.
lambpart = 0
for i in range(d):
    lambpart = lambpart + lamb*(gradu[i][i])
lambpart = sympy.simplify(lambpart, rational=True)
# lambpart = 0

stress = []
for i in range(d):
    scomp = []
    for j in range(d):
        s = 2*mu*0.5*(gradu[i][j] + gradu[j][i])
        if i == j:
            s = s + lambpart
        scomp.append(s)
    stress.append(scomp)

# We compute the divergence
div = []
for j in range(d):
    divcomp = 0
    for i in range(d):
        divcomp = divcomp + sympy.diff(stress[i][j], dervar[i])
    divcomp = sympy.simplify(divcomp, rational=True)
    div.append(divcomp)
