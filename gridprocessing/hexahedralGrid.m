function G = hexahedralGrid(P, H)
%Construct valid grid definition from points and list of hexahedra
%
% SYNOPSIS:
%   G = hexahedralGrid(P, T)
%
% PARAMETERS:
%   P     - Node coordinates.  Must be an m-by-3 matrix, one row for each
%           node/point.
%
%   T     - List of hexahedral corner nodes: an n-by-8 matrix where each
%           row holds node numbers for a hexahedron, with the following
%           sequence of nodes:
%                        (imax, jmin, kmin)
%                        (imax, jmax, kmin)
%                        (imax, jmax, kmax)
%                        (imax, jmin, kmax)
%                        (imin, jmin, kmin)
%                        (imin, jmax, kmin)
%                        (imin, jmax, kmax)
%                        (imin, jmin, kmax)
%
%
% RETURNS:
%   G     - Valid grid definition.
%
% EXAMPLE:
%
%   H = [1  2  3  4  5  6  7  8; ...
%        2  9 10  3  6 11 12  7; ...
%        9 13 14 10 11 15 16 12; ...
%       13 17 18 14 15 19 20 16; ...
%       17 21 22 18 19 23 24 20];
%
%   P = [1.2000       0  0.1860; ...
%        1.2000  0.0200  0.1852; ...
%        1.2000  0.0200  0.1926; ...
%        1.2000       0  0.1930; ...
%        1.1846       0  0.1854; ...
%        1.1844  0.0200  0.1846; ...
%        1.1848  0.0200  0.1923; ...
%        1.1849       0  0.1926; ...
%        1.2000  0.0400  0.1844; ...
%        1.2000  0.0400  0.1922; ...
%        1.1843  0.0400  0.1837; ...
%        1.1847  0.0400  0.1919; ...
%        1.2000  0.0601  0.1836; ...
%        1.2000  0.0600  0.1918; ...
%        1.1842  0.0601  0.1829; ...
%        1.1846  0.0600  0.1915; ...
%        1.2000  0.0801  0.1828; ...
%        1.2000  0.0800  0.1914; ...
%        1.1840  0.0801  0.1821; ...
%        1.1845  0.0800  0.1912; ...
%        1.2000  0.1001  0.1820; ...
%        1.2000  0.1001  0.1910; ...
%        1.1839  0.1001  0.1814; ...
%        1.1844  0.1001  0.1908];
%
%  G = hexahedralGrid(P, H);
%  G = computeGeometry(G);
%
%
% SEE ALSO:
%   `delaunay`, `tetrahedralGrid`, `grid_structure`

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


% Written by Jostein Natvig, SINTEF Applied Mathematics.

% NOTE:
% Implement unique face orientation given by magnitude of node indices:
% all faces have intrinsic orientation suth that the first two nodes in a
% face are the two smalles node indices of the face gemoetry.  This is
% needed to obtain a uniqe face representation before we find the
% face-to-cell map (G.faces.neighbors).

    % Definition of faces on hex
    Bot = H(:, [5, 6, 2, 1]);
    Top = H(:, [4, 3, 7, 8]);
    Est = H(:, [1, 2, 3, 4]);
    Wst = H(:, [8, 7, 6, 5]);
    Sth = H(:, [4, 8, 5, 1]);
    Nth = H(:, [2, 6, 7, 3]);

    % Rearrange such that for each row, the smalles node index is in the
    % first column
    nc       = size(H, 1);
    hfaces   = [Wst; Est; Sth; Nth; Bot; Top];
    [m, i]   = min(hfaces, [], 2); %#ok
    I        = bsxfun(@mod, bsxfun(@plus, (1:4)-1, i-1), 4)+1;
    k        = sub2ind(size(hfaces), repmat((1:nc*6)', [1, 4]), I);
    hfaces   = hfaces(k);


    % Rearrange such that for each row, the second smallest node is in the
    % second column.  We use 'i' as the 'sign' of each half-face.
    [m, i] = min(hfaces(:, [2,end]), [], 2); %#ok
    hfaces(i==2, :) = hfaces(i==2, [1,4,3,2]);

    % Better safe than sorry
    [m, k] = min(hfaces(:,[2,end]), [], 2); %#ok
    assert(all(k==1));
    clear k

    % id is [cellnumber, half-face tag]
    id       = [(1:nc)', repmat(1, [nc, 1]);...
                (1:nc)', repmat(2, [nc, 1]);...
                (1:nc)', repmat(3, [nc, 1]);...
                (1:nc)', repmat(4, [nc, 1]);...
                (1:nc)', repmat(5, [nc, 1]);...
                (1:nc)', repmat(6, [nc, 1]) ];

    % Sort rows to find pairs of cells sharing a face
    [hfaces, j] = sortrows(hfaces);

    % Run length compression to obtain unique face numbering
    [fnodes,n]  = rlencode(hfaces);
    nf          = numel(n);

    % Expand face numbers to each half-face
    N           = rldecode(1:nf,n,2)';

    % Accumulate results in face-2-cell map
    neigh  = accumarray([N, i(j)], id(j,1), [nf, 2]);

    % ... and a map of local directions
    hftag  = accumarray([N, i(j)], id(j,2), [nf, 2]);

    % ... and explicitly give number of nodes for each face.  So far, this
    % code is oblivious to pinch which may render faces triangular or
    % collapsed.
    nnodes = repmat(4, [nf, 1]);

    % Initialize grid structure
    G = emptyGrid(nc, P);

    % Add each face, given in terms of its nodes, cell neighbors.
    G = addFaces(G, reshape(fnodes', [], 1), nnodes,  neigh, hftag);

    % Record grid constructor in grid.
    G.type    = { mfilename };
    G.griddim = 3;
