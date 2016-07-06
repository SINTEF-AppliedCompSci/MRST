function [cg neighbors] = addVirtualCG(cg, neighbors, vmod)
% Add virtual cells to a coarse grid
%
% SYNOPSIS:
%   addVirtualCG(cg, neighbors)
%
% PARAMETERS:
%   cg        - coarse grid as defined by generateCoarseGrid(G, p)
%   neighbors - a cell array where cell i contains indices of the coarse
%               neighbors of coarse block i
%
%  RETURNS:
%   cg        - coarse grid with additional fields:
%                   - is_virtual containing bools indicating if a cell is
%                     virtual or not.
%                   - updated cells.num and cell.centroids for the new
%                   virtual blocks. The virtual blocks are not fully
%                   defined.
%   neighbors - updated neighbor structure, where virtual blocks are
%               neighbors to their parent (some coarse block on the edge)
%               as well as the virtual nodes grown from the neighbors of
%               the parent.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

    [edge_c ~] = find(cg.faces.neighbors == 0);
    Nc = numel(edge_c);

    Nf = cg.faces.num;
    cat_edge = zeros(Nf,1);
    for i = 1:Nf
        tmp = cg.cells.faces(:,1) == i;
        cats = cg.cells.faces(tmp,2);
        cat_edge(i) = max(cats);
    end
    newcents = zeros(Nc, 3);
    parents = zeros(Nc, 1);
    category = zeros(Nc,1);

    for i = 1:Nc
        % edge block this virtual blocks "grows" from
        eblock = max(cg.faces.neighbors(edge_c(i),:));
        % face centroid
        f_cent =     cg.faces.centroids(edge_c(i),:);
        % block centroid
        b_cent =     cg.cells.centroids(eblock   ,:);
        % add new centroids outside of the domain
        newcents(i,:) = b_cent + vmod*(f_cent - b_cent);
        % Set parent equal to the coarse block
        parents(i) = eblock;
        category(i) = cat_edge(edge_c(i));
    end

    for i = 1:Nc
        parents_neighbors = neighbors{parents(i)};
        globalindex = cg.cells.num+i;
        neigh = cg.cells.num + find(ismember(parents, parents_neighbors));
        % Remove self
        neigh = neigh(neigh~=globalindex);
        % If virtual neighbor is the same category, or is a direct neighbor
        % to the parent, add it
        neigh = neigh(category(neigh-cg.cells.num) == category(i) | parents(neigh-cg.cells.num) == parents(i));
        neighbors{globalindex} = [parents(i); neigh; globalindex]';
        neighbors{parents(i)} = [neighbors{parents(i)} globalindex];
    end

    cg.is_virtual = false(cg.cells.num + Nc, 1);
    cg.is_virtual(cg.cells.num+1 : end) = true;
    cg.cells.centroids = [cg.cells.centroids; newcents];
    cg.cells.num = cg.cells.num + Nc;


end
