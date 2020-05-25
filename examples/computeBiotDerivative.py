# compute external force for given displacement and pressure function

import sympy
from sympy import diff
from sympy import sin
from sympy import cos
from sympy import pi

sympy.init_printing(pretty_print=False)

x, y, z = sympy.symbols('x y z')
mu, lamb, alpha, K = sympy.symbols('mu0 lamb0 alpha K')
tau, rho = sympy.symbols('tau rho')

dervar = [x, y, z]

d = 2
if d == 3:
    u1 = y*(1 - x)*sin(2*pi*x*y)
    u2 = x*y**2*cos(x)
    u3 = x*y*z
    uvec = list([u1, u2, u3])
    p = u1
elif d == 2:
    u1 = y*(1 - x)*sin(2*pi*x*y)
    u2 = y**2*cos(x)
    uvec = list([u1, u2])
    p = u1

gradu = []

# gradu[i][j] gets the derivative of the i-th component with
# respect to the j-th variable
for i in range(d):
    g = []
    for j in range(d):
        g.append(diff(uvec[i], dervar[j]))
    gradu.append(g)

gradp = []
for i in range(d):
    gradp.append(diff(p, dervar[i]))

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

# We compute the divergence of sigma
divsigma = []
for j in range(d):
    divcomp = 0
    for i in range(d):
        divcomp = divcomp + diff(stress[i][j], dervar[i])
    divcomp = sympy.simplify(divcomp, rational=True)
    divsigma.append(divcomp)

# We compute the equation of conservation of momentum
momcons = []
for i in range(d):
    eq = divsigma[i] + alpha*gradp[i]
    eq = sympy.simplify(eq)
    momcons.append(eq)

# We compute the equation of mass conservation
divgradp = 0
for i in range(d):
    divgradp = divgradp + tau*K*diff(gradp[i], dervar[i])
divu = 0
for i in range(d):
    divu = divu + diff(uvec[i], dervar[i])
masscons = rho*p + alpha*divu + divgradp
masscons = sympy.simplify(masscons)
