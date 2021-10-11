function [qf, qf_vol] = calculateQF(G)
%
%
% SYNOPSIS:
%   function [qf, qf_vol] = calculateQF(G)
%
% DESCRIPTION:  Calculate elementary integrals that are used to assemble the
% stiffness matrix for the 2D case. 
%
% PARAMETERS:
%   G - Grid structure
%
% RETURNS:
%
%   qf     -    Elementary assembly integrals : One (2D) vector value in each cell,
%               which corresponds to the two components of the integral of the
%               basis function in each coordinate over the cell faces (see (74) in
%               [Gain et al: doi:10.1016/j.cma.2014.05.005], then faces there
%               correspond to edges here. We do not divide by 2*(cell volume) here).
%   qf_vol -    Elementary assembly integrals : one scalar value for each
%               node, wich corresponds to the weights that are used to
%               compute th L^2 projection, see VEM_linElast.m
%
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


    assert(G.griddim == 2); 
    
    qf = zeros(size(G.cells.nodes, 1), 2); 
    cellno = rldecode([1 : G.cells.num]', diff(G.cells.nodePos)); % #ok

    cells   = 1 : G.cells.num;
    lcells = rldecode(cells', diff(G.cells.nodePos)');
    
    % For each cell, indices of the first to the second-to-last node 
    inodes1 = mcolon(G.cells.nodePos(cells), G.cells.nodePos(cells + 1) - 2)'; 
    
    % For each cell, indices of the second to the last node 
    inodes2 = mcolon(G.cells.nodePos(cells) + 1, G.cells.nodePos(cells + 1) - 1)'; 

    % For each cell, indices of each face 
    ifaces  = mcolon(G.cells.facePos(cells), G.cells.facePos(cells + 1) - 1)'; 
    
    % For each cell, indices of each face 
    faces   = G.cells.faces(ifaces, 1); 
    
    % Orienting normals ('N') so that they always point out of the current cell and
    % into the neighbor cell
    sign    = 2 * (G.faces.neighbors(faces, 1) == cellno) - 1; 
    N       = bsxfun(@times, G.faces.normals(faces', :), sign); 

    % For each cell node, add up the (scaled) normals of the two adjacent faces and
    % divide by two.
    relvec = G.faces.centroids(faces, :) - G.cells.centroids(lcells, :);
    tetvols = sum(N.*relvec, 2);
    qf_vol = zeros(numel(G.cells.nodes), 1);
    qf_vol(inodes1) = qf_vol(inodes1, :) + tetvols(inodes1); 
    qf_vol(inodes2) = qf_vol(inodes2) +  tetvols(inodes1);
    qf_vol(G.cells.nodePos(cells)) = qf_vol(G.cells.nodePos(cells)) + ...
        tetvols(G.cells.nodePos(cells + 1) - 1);
    qf_vol(G.cells.nodePos(cells + 1) - 1) = qf_vol(G.cells.nodePos(cells + 1) - 1) + ...
        tetvols(G.cells.nodePos(cells + 1) - 1);                                    
    qf_vol = qf_vol/4;
    
    qf(inodes1, :) = qf(inodes1, :) + N(inodes1, :); 
    qf(inodes2, :) = qf(inodes2, :) + N(inodes1, :);
    qf(G.cells.nodePos(cells), :) = qf(G.cells.nodePos(cells), :) + ...
        N(G.cells.nodePos(cells + 1) - 1, :); 
    qf(G.cells.nodePos(cells + 1) - 1, :) = qf(G.cells.nodePos(cells + 1) - 1, :) + ...
        N(G.cells.nodePos(cells + 1) - 1, :);
    qf = qf / 2;

end

