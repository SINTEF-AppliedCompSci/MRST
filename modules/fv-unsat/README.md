# fv-unsat: An MRST module for unsaturated poroelasticity

## Description
This module implements a discretization of the equations of unsaturated poroelasticity using cell-centered finite volume methods, specifically the Multipoint Flux Approximations (MPFA) for flow and the Multipoint Stress Approximation (MPSA) for (poro)-elasticity. If mechanical effects are neglected, the set of equations reduces to the well-known Richards' equation. The module is written based on the Automatic Differentiation framework provided by MRST.

There are four numerical examples included in this module:
* convAnalysisRE.m
* convAnalysisUnsatBiot.m
* waterInfiltrationRE.m
* desiccationUnsatBiot.m

The first two are convergence tests and the last two are practical applications. Even though the numerical tests are well documented, they are not meant as tutorials, but rather included for demonstrative purposes. To learn the basics regarding the module usage, we recommend starting with `waterInfiltrationRE.m`.

This module was largely based on the Master's thesis:
* Varela, Jhabriel. Implementation of an MPFA/MPSA-FV Solver for the Unsaturated Flow in Deformable Porous Media. MS thesis. The University of Bergen, 2018.

For an introduction to MPFA, refer to:
* Aavatsmark, Ivar. "An introduction to multipoint flux approximations for quadrilateral grids." Computational Geosciences 6.3-4 (2002): 405-432.

For an introduction to MPSA, refer to:
* Keilegavlen, Eirik, and Jan Martin Nordbotten. "Finite volume methods for elasticity with weak symmetry." International Journal for Numerical Methods in Engineering 112.8 (2017): 939-962.
* Nordbotten, Jan Martin. "Cell‐centered finite volume discretizations for deformable porous media." International journal for numerical methods in engineering 100.6 (2014): 399-418.
* Nordbotten, Jan Martin. "Stable cell-centered finite volume discretization for Biot equations." SIAM Journal on Numerical Analysis 54.2 (2016): 942-968.

You can also check:
* Nordbotten, J. M., & Keilegavlen, E. (2021). An introduction to multi-point flux (MPFA) and stress (MPSA) finite volume methods for thermo-poroelasticity. In Polyhedral methods in geosciences (pp. 119-158). Cham: Springer International Publishing.

## Requirements
* MRST (Tested version: 2019b)
* MATLAB (Tested version: R2019a)

## MRST dependencies
* [fvbiot](https://github.com/pmgbergen/fvbiot) (Please see Troubleshooting section below)
* [upr](https://www.sintef.no/projectweb/mrst/modules/upr/)

## Troubleshooting
* If you are using an MRST version <= 2019b, you will have to clone the [fvbiot](https://github.com/pmgbergen/fvbiot) repository, and replace it manually in your MRST folder. Note that fvbiot is located inside the "modules" folder.

## Cite
If you use `fv-unsat`, please cite:

Varela, J., Gasda, S. E., Keilegavlen, E., & Nordbotten, J. M. (2021). A finite-volume-based module for unsaturated poroelasticity. Advanced Modeling with the MATLAB Reservoir Simulation Toolbox.

## Contact
Jhabriel Varela (jvarela@uptp.edu.py).

