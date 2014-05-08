function rock=setupRockStructure(G, perm, poro)
% Sets up a rock structure for grid G, with uniform permeability
% and porosity as specified by 'perm' and 'poro'.
% 
% SYNOPSIS:
%   rock=setupRockStructure(G, perm, poro)
%
% PARAMETERS:
%   G     - the grid for which the rock structure will be generated
%   perm  - uniform permeability to be assigned
%   poro  - uniform porosity to be assigned
%
% RETURNS:
%   rock  - the generated rock structure.

num_cells = G.cells.num;

rock.perm = repmat(perm, num_cells, 1);
rock.poro = repmat(poro, num_cells, 1);

