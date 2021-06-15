function out = getMexDiscreteDivergenceJacPrecomputes(model)
%Compute Mapping Primitives for Accelerated Discrete Divergence Operator
%
% SYNOPSIS:
%   map = getMexDiscreteDivergenceJacPrecomputes(model)
%
% PARAMETERS:
%   model - Discretised physical model including neighbourship information
%
% RETURNS:
%   map - Mapping primitives.  Structure with the following fields:
%           - facePos   - (#cells + 1)-by-1 indirection map into connection
%                         information.  Plays the same role as the field
%                         with the same name in the `grid_structure`.
%
%           - faces     - Zero-based index of cells' connecting interfaces.
%                         The connecting faces of cell 'i' are
%
%                            faces(facePos(i) : facePos(i + 1) - 1)
%
%           - cells     - One-based index of connecting/neighbouring cells.
%                         The connecting/neighbouring cells of cell 'i' are
%
%                            cells(facePos(i) : facePos(i + 1) - 1)
%
%                         Moreover, the 'cells' array encodes the direction
%                         of flow through a sign convention.  In particular
%                         the neigbouring cell index is listed as POSITIVE
%                         if the neighbour is a source of flow INTO cell
%                         'i' and NEGATIVE if the flow is nominally OUT of
%                         cell 'i'.
%
%           - cellIndex - Zero-based non-zero position of diagnoal matrix
%                         element corresponding to each cell.  Integer array
%                         of size #cells-by-1.
%
% NOTE:
%   This function uses SORTROWS.
%
% SEE ALSO:
%   `mexDiscreteDivergenceJac`, `grid_structure`, sortrows.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    nc = model.G.cells.num;
    N  = model.operators.N;

    f = (1 : size(N, 1)) .';

    % [ self, other, faceID ]
    cell_connections = sortrows([[N ; fliplr(N)], [f ; f]]);
    sum_to_cell = @(x) accumarray(cell_connections(:, 1), x, [nc, 1]);

    % Positive if 'self' is in second column of N, negative otherwise.
    sign = 2*(N(cell_connections(:, 3), 1) ~= cell_connections(:, 1)) - 1;

    localCellIndex = ...
       sum_to_cell(cell_connections(:, 2) < cell_connections(:, 1));

    out = struct('facePos'  , cumsum([0 ; sum_to_cell(1)]) , ...
                 'faces'    , cell_connections(:,3)  - 1   , ...
                 'cells'    , sign .* cell_connections(:,2), ...
                 'cellIndex', localCellIndex);
end
