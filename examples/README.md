

# Convergence Tests

-   The directory [`convergencetests`](convergencetests) contains scripts for convergence tests.


# MPSA - Mechanics


## [`assemblyMpsaExample.m`](assemblyBiotExample.m)

-   Example using bare assembly functionalities.
-   We use directly the assembly output matrices to set up and solve the system.


## [`tiltedExample.m`](tiltedExample.m)

-   Example using sliding condition in a non-Cartesian direction.


## [`mpsaExample.m`](mpsaExample.m)

-   We use [`MechModel`](../models/MechModel.m) AD model to set up the equations.
-   We use predefined mechanical test cases ('2d-linear', '2d-refinement', etc.).


# MPFA - Flow


## [`assemblyMpfaExample.m`](assemblyMpfaExample.m)

-   Example using bare assembly functionalities.
-   With this example, we verify on a skewed grid that a linear pressure field is computed exactly by MPFA.


## [`assembleMpfaExample2.m`](assembleMpfaExample2.m)

-   Example using bare assembly functionalities
-   We use the grids from predefined mechanical test cases ('2d-linear', '2d-refinement', etc), as in `mpsaExample`


## [`mpfaBlackoilExample.m`](mpfaBlackoilExample.m)

-   Example that demonstrates grid effects for TPFA and how MPFA remove those.
-   We consider an oil-water system and inject water in a oil reservoir. The setup is such that the analytical solution is
    symmetric. The symmetry of the solution is clearly broken by TPFA.
-   We use [`MpfaBlackOilModel`](../models/MpfaBlackOilModel.m) AD model.


# Biot MPSA-MPFA - Poroelasticity


## [`assemblyBiotExample.m`](assemblyBiotExample.m)

-   Example using bare assembly functionnalities.
-   We use the same mechanical test cases ('2d-linear', '2d-refinement', etc) as in `mpsaExample`.


## [`mandel.m`](mandel.m)

-   Mandel test case.
-   We use a the AD model [`MandelModel`](file:///home/xavier/Matlab/Projects/project-mpsaw/models/MandelModel.m) which has been
    specifically implemented to deal with the special boundary conditions of the Mandel test case.


## [`topforceExample.m`](topforceExample.m)

-   Test case where a controlled force is applied at the top.
-   The intensity of the force, which varies in time, is provided through a standard MRST schedule.
-   We use [`BiotModel`](../models/BiotModel.m) AD model.


## [`biotBlackoilExample.m`](biotBlackoilExample.m)

-   Example where MPSA-MPFA/TPFA is coupled with the object-oriented automatic-differentiation (OO-AD) framework of MRST.
-   Same setup as example `mpfaBlackoilExample` showing how the methods cope with grid effects.
-   We use [`BiotBlackOilModel`](../models/BiotBlackOilModel.m) and [`BiotTpfaBlackOilModel`](../models/BiotTpfaBlackOilModel.m) AD model.


## [`biotCompositionalExample.m`](biotCompositionalExample.m)

-   A compositional simulation (two components: CO2 and n-decane, two phases: gas and liquid, cubic equation of state) which
    also couples MPSA-MPFA with the `compositional` module of MRST
-   We use [`BiotCompositionalOilModel`](../models/BiotCompositionalModel.m) and [`BiotTpfaCompositionalOilModel`](../models/BiotTpfaCompositionalModel.m) AD model.

