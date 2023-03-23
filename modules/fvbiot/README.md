This package implements a discretization of poro-mechanics by cell centered 
finite volume methods. 

The package provides discretization of three different equations:
   1. Scalar elliptic equations (Darcy flow), using Multipoint flux approximations. 
   2. Linear elasticity, using Multipoint stress approximations
   3. Poro-mechanics, e.g. coupling terms for the combined system of 1. and 2.

Examples of usage for the three options can be found in mpfa_ex, mpsa_ex and biot_ex, respectively.

Further information on the methods implemented herein, in particular the discretization of elasticity and poro-mechanics, can be found in the following papers (all open access):

J. M. Nordbotten: Stable cell-center finite volume disrcetizaiton for
   Biot equations - SIAM J. Numer. Anal. 54(2) 942-968.

J. M. Nordbotten: Convergence of a cell-centered finite volume discretization 
   for linear elasticity - SIAM J. Numer. Anal. 53(2) 2605-2625.

E. Keilegavlen, J. M. Nordbotten: Finite volume methods for elasticity with weak
   symmetry - Int. J. Num. Meth. Eng, 2017, doi: 10.1002/nme.5538.
   
The implementation of elasticity is based on the forumlation described in the latter. 
For information on the MPFA discretizations for the flow equation, see 
   I. Aavatsmark: An introduction to multipoint flux approximations for quadrilateral grids, Comput. Geosci. 6(3-4) 405-432.

The is compatible with the (freely available) Matlab Reservoir Simulation
Toolbox (MRST) provided by SINTEF ICT, see http://www.sintef.no/projectweb/mrst/
In particular, the grid data structure is that of MRST, and it is assumed that 
MRST is in the matlab path. 

The code has been tested with Matlab R2015b, and mrst-2016a.
