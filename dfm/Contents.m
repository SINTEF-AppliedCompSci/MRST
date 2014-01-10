% Discrete Fracture Matrix (DFM) module
%
% Routines supporting the DFM method. Examples in /examples
% Add module by typing MrstModule add dfm
% AUTHOR: tor.harald.sandve@iris.no
%
% new files
%   addhybrid              - Generates the hybrid cells
%   computeHybridTrans     - Computes the hybrid transmissibilities using TPFA
%   computeMultiPointTrans - Computes the hybrid transmissibilities using
%                            O-MPFA
%   testNormals            - Rutine to check that the direction of the normals
%                            is from neighbor 1 to 2.
%   /plotting/plotFractures - Plots 2d hybridcells.
%   build_fractures_mod     - read fracture date from open-office or .mat files
%   plotEdges 	 		- plot lines
%   nodeType			- return the node type: 0: interior, 1: edge 2: corner.
%   triangulate			- create Delaunay triangulation based on a cloud of points and a set of constrains
%
% Files modified from  core MRST functions. Use these to
%
%   computeTrans_DFM           - Modified to compute hybrid-normal
%                                transmissibilities using a modified TPFA
%   computeMultiPointTrans_DFM - Modified to compute hybrid-normal
%                                transmissibilities using a modified MPFA
%   computeTPFA_DFM            - Modified to allow for cell2cell connections
%   computeMPFA_DFM            - Modified to allow for cell2cell connections
%   twophaseJacobian_DFM       - Modified to allow for cell2cell connections
%   explicitTransport_DFM      - Modified to allow for cell2cell connections
%   implicitTransport_DFM      - Modified to allow for cell2cell connections
%   removeInternalBoundary_DFM - Modified to update grid geometry as well as the topology
%   /plotting/                 - Modification to plotting rutines to handle
%  computeTimeOfFlight		- Modified to allow for cell2cell connections
%
%
% /private.
%  - Files needed since we can not access private folders from this location
%
% /examples:
%  matlabFractureGrid   - Unstructured fracture grid created by Matlab functions
%  triangleFractureGrid - Unstructured fracture grid created using the Triangle software
%
