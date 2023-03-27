function  G = readIGEMSSurface(location, name, i, z_levels, z_height, coarse)
% Read one of the IGEMS surfaces from its IRAP file.
% 
% SYNOPSIS:
%   G = readIGEMSSurface(location, name, i, z_levels, z_height, coarse)
%
% PARAMETERS:
%   location - location of the top directory containing the individual IGEMS
%              surface subdirectories (flatxxx, FMMxxx, OSSxxx, etc.)
%   name     - name of the surface type (flatNP1, FMMUP2, OSS, etc.)
%   i        - surface realization number (an integer between 1 and 99
%              inclusive)
%   z_levels - number of levels to establish in the z-direction
%   z_height - total height in z-direction (in meters)
%   coarse   - 2-component vector of integers.  Default is [1 1].  If
%              different from default, specifies a downsampling factor in the
%              number of x-cells and y-cells in the grid.
%
% RETURNS:
%   G        - the IGEMS surface as a 3D MRST-grid 

% @@ Requires the mex/libgeometry module 

xy_cellsz = 100; % expected cellsize, in meters, for x and y direction
if (isempty(coarse)) coarse = [1 1]; end % no coarsening by default

% Reading grid information from file
[x, y, Z, angle] = readIrapClassicAsciiSurf([location, '/', name, '/', ...
                                             int2str(i), '_', name, '.irap']); %#ok
% All horizontal increments should be 100 m
assert(all(diff(x) == xy_cellsz));
assert(all(diff(y) == xy_cellsz));

% Downsampling, if necessary
Z = Z(1:coarse(1):end, 1:coarse(2):end);  % unchanged when coarse=[1 1]
num_cells_xy = size(Z) - 1;
phys_size = [num_cells_xy, 1] .* [coarse, 1] .* [xy_cellsz xy_cellsz z_height];

% Generating cartesian 3D-grid, and add Z to the node heights
G = cartGrid([num_cells_xy, z_levels], phys_size);
G.nodes.coords(:,3) = G.nodes.coords(:,3) + repmat(Z(:), z_levels+1, 1);

% Computing the midpoint of the top surface of each cell (average of corners)
Z_cells = 1/4 * (Z(1:end-1, 1:end-1) + Z(2:end, 1:end-1) + ...
                 Z(1:end-1, 2:end)   + Z(2:end, 2:end));

% Removing invalid cells, and computing resulting geometry
broken_cells = repmat(isnan(Z_cells), [1 1 z_levels]);
%G = mcomputeGeometry(removeCells(G, find(broken_cells)));
G = computeGeometry(removeCells(G, find(broken_cells)));









