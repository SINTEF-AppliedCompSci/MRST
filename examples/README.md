

# Convergence Tests

-   The directory [`convergencetests`](convergencetests.m) contains all the convergence tests.


# MPSA - Mechanics


## [`assemblyMpsaExample`](assemblyBiotExample.m)

-   Example using bare assembly functionalities.


## [`tiltedExample`](tiltedExample.m)

-   Example using sliding condition in a non-Cartesian direction.


## [`mpsaExample`](mpsaExample.m)

-   use [`MechModel`](../models/MechModel.m) AD model.
-   use predefined mechanical test cases ('2d-linear', '2d-refinement', etc.).


# MPFA - Flow


## [`assemblyMpfaExample`](assemblyMpfaExample.m)

-   Example using bare assembly functionalities.


## [`assembleMpfaExample2`](assembleMpfaExample2.m)

-   Example using bare assembly functionalities
-   use grids from predefined mechanical test cases ('2d-linear', '2d-refinement', etc)


## [`mpfaBlackoilExample`](mpfaBlackoilExample.m)

-   Example that demonstrates grid effects for TPFA and how MPFA remove those.
-   use [`MpfaBlackOilModel`](../models/MpfaBlackOilModel.m) AD model.


# Biot MPSA-MPFA - Poroelasticity


## [`assemblyBiotExample`](assemblyBiotExample.m)

-   use bare assembly functionnalities.
-   uses mech test cases ('2d-linear', '2d-refinement', etc).


## [`mandelExample`](mandelExample.m)

-   Mandel test case.
-   use [`BiotModel`](../models/BiotModel.m) AD model.


## [`biotBlackoilExample`](biotBlackoilExample.m)

-   Same setup as example `tiltedExample` showing how method copes with grid effects
-   use [`BiotBlackOilModel`](../models/BiotBlackOilModel.m) AD model.


## [`biotExample`](biotExample.m)

-   same setup as mandel test case but using wells and a fixed external force at the top
-   use [`BiotModel`](../models/BiotModel.m) AD model

