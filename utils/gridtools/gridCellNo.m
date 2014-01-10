function cellno = gridCellNo(G, varargin)
% Find map from half-faces to cells given a set of cells and a grid G
%
% SYNOPSIS:
%   cellno = gridCellNo(G)
%   cellno = gridCellNo(G, c)
%
% PARAMETERS:
%   G    - Grid structure
%   c    - Cells for which the faces are required
%
%
% RETURNS:
%   cellno - One element for each half-face. cellno(hf) gives the cell
%            where halfface hf belongs.
%

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


    if nargin == 1
        cells = 1:G.cells.num;
    else
        cells = varargin{1};
    end

    % Find number of faces per cell
    nf = diff([G.cells.facePos(cells), G.cells.facePos(cells+1)], [],2);

    % Find the cell index of each face
    cellno = rldecode(1:numel(cells), nf, 2) .';
end
