% Multiscale Finite Volume Method for Discrete Fracture Matrix models (MSFV-DFM) module
% The module depends on the MSFVM and the DFM module
%
% Routines supporting the MSFVM-DFM module. Examples in /examples
% Add module by typing MrstModule add MSFV-DFM
% AUTHOR tor.harald.sandve@iris.no
%
% new files
%   createGridHierarchy    	- Generates the fine-scale grid and the dual-coarse grid
%
%   /util/findClosestCell 	- Find the closest cells to a given coordinate in 2D
%   /util/getFractureCells 	- Extract cell indexes corresponding to fractures and fracture intersection
%				with a given tag
%   /util/makeCartFrac 		- Create a net of Cartesian fractures
%   /util/mergeCells    	- merge two cells that shares a face
%   /util/mergeFaces    	- merge two faces that have the same neighbors
%   /util/mergeFractures 	- merge two fracture sets
%   /util/removeDuplicates 	- remove duplicated fractures and vertices
%   /private/makeCoarseDualGrid - Create the Coarse dual grid
%   /private/makeFineGridHybrid - Create a fine-grid constrained by fractures and the dual coarse grid
%				  The constrains are represented by hybrid cells of first kind
%				  and their intersection with hybrid cells of second kind
%   /private/readFractures   	- Extract fracture network from files
%
% Files modified from the MSFVM module.
%
%   solveMSFV_TPFA_Incomp_DFM  	- Modified to allow for hybrid cells
%   /util/plotDual_DFM   	- Modified to plot hybrid cells
%   /util/removeCells_mod  	- Modified to update additional face fields as tags, centroids etc.
%
%   partitionDualDFM          - Create the dual grid structure as used by MSFVM
%   partitioningByAggregation - Create primal partitioning by aggregation with node cells as seeds.
