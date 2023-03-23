

# Implementation of MPSA with weak symmetry for linear elasticity and poroelasticity

This module provides an implementation of the Multiple Point Stress Approximation (MPSA) using weak symmetry. The method is also
coupled with MPFA to solve poroelasticity problems.


## Reference papers

The implementation follows the description done in the two following papers:

-   *Finite volume methods for elasticity with weak symmetry* Keilegavlen, Eirik and Nordbotten, Jan Martin, International Journal for Numerical Methods in Engineering 2017 <https://doi.org/10.1002/nme.5538>
-   *Stable cell-centered finite volume discretization for Biot equations*, Nordbotten, Jan Martin, SIAM Journal on Numerical Analysis 2016, <https://doi.org/10.1137/15M1014280>


## Installation

-   The MPSAW module is included in every MRST release as a *snapshot* of this repository at the time of the
    release. Therefore, the easiest way to install the module is to install the latest release of MRST which can be
    found [here](https://www.sintef.no/projectweb/mrst/download/) with the installation instructions.
-   For developpers, it is also possible to use directly the module from this repository. MRST should be installed
    (either a release version or a [development version](https://bitbucket.org/mrst/mrst-core/wiki/Home)). Then, the module should be registered as described in an
    example [here](https://bitbucket.org/mrst/mrst-core/wiki/Home) (end of page).


## Examples

Examples can be found in [`examples`](examples/) directory

