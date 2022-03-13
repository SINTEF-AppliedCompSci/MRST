function cell_field = nodeFieldToCellField(Gt, node_field)
%
% Convert a field of node-based values to a field of cell-based values.
% Each node gets the value associed with the node with lowest z-value of its four
% corner nodes.
% This function was primarily written for projecting the spill-regions computed 
% based on nodes, onto cells.
% 
% SYNOPSIS:
% cell_field = nodeFieldToCellField(Gt, node_field) 
%
% PARAMETERS:
% Gt         - 2D grid for which the field of node-based values has been defined.
% node_field - the node-based field of values (e.g. spill region number)
% 
% RETURNS:
% cell_field - the cell-based field resulting from the projection of 'node_field' down
%              to the associated grid cells.
% 
% SEE ALSO: 
% `computeNodeTraps`

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    % Determining the corner of each cell with lowest depth (i.e. situated 'highest')
    cnodes = activeCellNodes(Gt);
    num_cells = size(cnodes, 2);
    [min_z, min_z_ix ] = min(Gt.nodes.z(cnodes)); %#ok backwards compatability
    
    % Choosing the node_field value of the chosen (i.e. 'highest') corner of each cell 
    % to represent the value for that cell
    cell_field = node_field(cnodes(sub2ind(size(cnodes),min_z_ix, 1:num_cells)));
end
