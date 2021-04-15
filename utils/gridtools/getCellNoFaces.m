function [cellNo, cellFaces, isNNC] = getCellNoFaces(G)
%Get a list over all half faces, accounting for possible NNC
%
% SYNOPSIS:
%   [cellNo, cellFaces, isNNC] = getCellNoFaces(G)
%
% DESCRIPTION:
%   This utility function is used to produce a listing of all half faces
%   in a grid along with the respective cells they belong to. While
%   relatively trivial for most grids, this function specifically accounts
%   for non-neighboring connections / NNC.
%
% REQUIRED PARAMETERS:
%   G         - Grid structure with optional .nnc.cells field.
%
% RETURNS:
%   cellNo    - A list of length number of geometric half-faces + 2* no. nnc
%             connections where each entry corresponds to the cell index of
%             that half face.
%
%   cellFaces - A list of length number of geometric half-faces + 2* number
%             of nnc connections where each entry is the connection index.
%             For the first entries, this is simply the face number.
%             Otherwise, it is the entry of the NNC connection.
%
%   isNNC     - A list with the same length as cellNo / cellFaces,
%               containing a boolean indicating if that specific connection
%               originates from a geometric face or a NNC connection
% SEE ALSO:
%   `rldecode`

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

    % Find the cell index of each face for each cell
    cellNo   = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
    % List of faces per cell with same size as cellNo
    cellFaces = G.cells.faces(:,1);
    
    % Mapping to show which entries in cellNo/cellFaces are resulting from
    % NNC and not actual geometric faces.
    isNNC = false([numel(cellNo), 1]);
    % If NNC is present, we add these extra connections to cellno &
    % cellFaces to allow consistent treatment.
    if isfield(G, 'nnc') && isfield(G.nnc, 'cells')
        nnc_cells = [G.nnc.cells(:, 1); G.nnc.cells(:, 2);];
        % NNC is located at the end after the regular faces
        nnc_faceno = G.faces.num + (1:size(G.nnc.cells, 1)).';
        cellNo = [cellNo; nnc_cells];
        cellFaces = [cellFaces; nnc_faceno; nnc_faceno];
        % Added connections are NNC
        isNNC = [ isNNC ; true([numel(nnc_cells), 1]) ];
    end
end