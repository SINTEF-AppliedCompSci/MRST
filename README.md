This package implements a discretization of the (three-dimensional) incompressible Richards' equation by 
cell centered finite volume methods, specifically Multi-Point Flux Approximation. 

Due to the non-linear nature of the water retention curves (van Genuchten-Mualem), the resulting set of 
partial differential equations is also non-linear. To approximate the solutions at each time step, 
we use a classical Newton iterative scheme, where the Automatic Differentiation capabilities
of MRST are exploited. 

Examples of usage for three scenarios can be found in waterInfiltration1D, waterInfiltration3D and waterEvaporation3D, respectively.

For information on the MPFA discretizations for the flow equation, see 
   I. Aavatsmark: An introduction to multipoint flux approximations for quadrilateral grids, Comput. Geosci. 6(3-4) 405-432.

The codes are compatible with the (freely available) Matlab Reservoir Simulation
Toolbox (MRST) provided by SINTEF ICT, see http://www.sintef.no/projectweb/mrst/
In particular, the grid data structure is that of MRST, and it is assumed that 
MRST is in the MATLAB path.

Also, it is assumed that the fvbiot package (freely available) which provides the discretization of the flux is in the MATLAB path as well. The fvbiot package was developed by Eirik Keilevgalen, Ph.D. from the Porous Media Group of the University of Bergen and can be downloaded from https://github.com/pmgbergen/fvbiot

The code has been tested with Matlab R2017b, and mrst-2016a.

We highly recommend you to use the MATLAB publish tool, where you will be able to go through each section of the scripts, which were carefully documented.

If you use RE-MPFA, please cite:
Varela, J. (2018). Implementation of an MPFA/MPSA-FV Solver for the Unsaturated Flow in Deformable Porous Media (Master's thesis, The University of Bergen).

You can report bugs or make consultations by e-mailing: 
Jhabriel Varela, M.Sc. (jhabriel@gmail.com).
