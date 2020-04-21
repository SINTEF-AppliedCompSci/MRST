# fv-unsat: An MRST module for unsaturated poroelasticity

## Description
This module implements a discretization of the (three-dimensional) equations of unsaturated poroelasticity using
cell-centered finite volume methods, specifically MPFA-O and MPSA-W. Note that if  mechanical effects are neglected, the set of equations reduce to the well-known Richards' equation. The module is written based on the Automatic Differentiation framework provided by MRST.

There are four numerical examples accompanying this module:
* convAnalysisRE.m
* convAnalysisUnsatBiot.m
* waterInfiltrationRE.m
* desiccationUnsatBiot.m

The first two are convergence tests and the last two are practical applications. Even though the numerical tests are well documented, they are not meant as tutorials, but rather included for demonstrative purposes. To learn the basics regarding the module usage, we recommend waterInfiltrationRE.m.

This module was largely based on:
* Varela, Jhabriel. Implementation of an MPFA/MPSA-FV Solver for the Unsaturated Flow in Deformable Porous Media. MS thesis. The University of Bergen, 2018.

For an introduction to MPFA:
* Aavatsmark, Ivar. "An introduction to multipoint flux approximations for quadrilateral grids." Computational Geosciences 6.3-4 (2002): 405-432.

For an introduction to MPSA:
* Keilegavlen, Eirik, and Jan Martin Nordbotten. "Finite volume methods for elasticity with weak symmetry." International Journal for Numerical Methods in Engineering 112.8 (2017): 939-962.
* Nordbotten, Jan Martin. "Cell‐centered finite volume discretizations for deformable porous media." International journal for numerical methods in engineering 100.6 (2014): 399-418.
* Nordbotten, Jan Martin. "Stable cell-centered finite volume discretization for Biot equations." SIAM Journal on Numerical Analysis 54.2 (2016): 942-968.
 
## Requirements
* MRST (Tested version: 2019b)
* MATLAB (Tested version: R2019a)

## MRST dependencies
* [fvbiot](https://github.com/pmgbergen/fvbiot) (Please see Troubleshooting section below)
* [distmesh](http://persson.berkeley.edu/distmesh/)

## Troubleshooting
* If you are using an MRST version <= 2019b, you will have to clone the [fvbiot](https://github.com/pmgbergen/fvbiot) repository, and replace it manually in your MRST folder. Note that fvbiot is located inside the "modules" folder.

## Cite
If you use fv-unsat, please cite:
* *CHAPTER GOES HERE*

## Contact
Jhabriel Varela (jhabriel.varela@uib.no).

