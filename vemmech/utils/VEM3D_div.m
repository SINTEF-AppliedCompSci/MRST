function [div] = VEM3D_div(G)
% Discrete divergence operator for the virtual element method in 3D
%
% SYNOPSIS:
%   function [div] = VEM3D_div(G)
%
% DESCRIPTION: Computes a discrete divergence operator in 3D. This discrete div
% operator is a mapping from node-valued displacement vector to cell-valued 2D
% vector. Node-valued displacement vectors correspond to the degrees of freedom
% that determine for each cell a displacement function over the cell via the
% virtual basis function. The discrete divergence operator that is assembled here
% computes the L^2 projection, cell-wise, of this displacement function. For
% more detail, see paper [Andersen et al: http://arxiv.org/abs/1606.09508v1].
%
% PARAMETERS:
%   G - Grid structure
%
% RETURNS:
%   div - matrix corresponding to the discrete divergence operator.
%
% EXAMPLE:
%
% SEE ALSO:
%

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


    assert(G.griddim == 3);
    cells     = 1:G.cells.num;
    ifaces    = mcolon(G.cells.facePos(cells), (G.cells.facePos(cells + 1) - 1));
    faces     = G.cells.faces(ifaces, 1);
    inodes    = mcolon(G.faces.nodePos(faces), (G.faces.nodePos(faces + 1) - 1));
    nodes     = G.faces.nodes(inodes);
    cellfaces = rldecode(cells', diff(G.cells.facePos));
    nvec      = bsxfun(@times, G.faces.normals(faces, :), ...
                       (2 * (G.faces.neighbors(faces, 1) == cellfaces) - 1));
    nfl       = G.faces.nodePos(faces + 1) - G.faces.nodePos(faces);
    nvec  = bsxfun(@rdivide, nvec, G.faces.areas(faces));
    [qc, qf] = calculateQC(G);
    matel = reshape(bsxfun(@times, rldecode(nvec, nfl), qf(inodes))', [], 1);
    jind  = mcolon(G.griddim * (nodes - 1) + 1, G.griddim * (nodes - 1) + G.griddim);
    iind  = rldecode(cellfaces, nfl * 3);
    div   = sparse( iind', jind', matel(:), G.cells.num, G.nodes.num .* G.griddim);

end
