function [upneigh, nhood] = findUpslopeNeighbor(xyz, neighs, regions)
%
% Determine cell (or node) neighbors, and determine steepest upslope
%
% SYNOPSIS:
%   function [upneigh, nhood] = findUpslopeNeighbor(xyz, neighbors, regions)
%
% DESCRIPTION:
%   Determine cell (or node) neighbors, and determine steepest upslope.
%   Neighbors include diagonal neighbors and the cell (or node) itself, so
%   the max. number of neighbors is 9.
% PARAMETERS:
%   xyz     - 3D coordinates of cells (or nodes)
%   neighs  - (n,2) list of the 'n' neighbor relations, where each line
%             gives the indices of two cells that are neighbors
%   regions - optional parameter that assigns a 'region' number to each
%             cell (or node).  Neighborhood relations across regions are
%             eliminated.
%
% RETURNS:
%   upneigh - most upslope neighbor for each cell (or node).  Can be the cell
%             itself in case of a sommet node.
%   nhood   - a matrix giving the complete neighborhood of each cell (or
%             Each line represents a node, and contains the indices of nodes
%             in its neighborhood.  Number of columns is governed by the
%             largest neighborhood.  If a cell has fewer neighbors than the
%             number of columns, the surplus entries are given the value zero.
% EXAMPLE:
%
% To compute the cell neighborhood and upslope neighbors of cells in a top
% surface grid 'Gt':
%   [upslope, neigh] = findUpslopeNeighbor([Gt.cells.centroids, Gt.cells.z],
%                                          Gt.faces.neighbors);
% 

    num = size(xyz,1); % number of cells (or nodes)
    
    % Remove neighbor relations representing exterior faces
    neighs = neighs(find(prod(neighs, 2)~=0),:);
    
    % Producing adjacency matrix for the 4-neighborhood of each cell
    adj = sparse(neighs(:,1), neighs(:,2), 1, num, num, 2*size(neighs,1));
    adj = adj + adj';

    % Include diagonal neighbors and cell itself
    [I,J] = ind2sub(size(adj), find(adj*adj==2));
    adj = spones(adj + sparse(I, J, 1, num, num) + speye(num));
    
    % Producing list of cell neighbor pairs
    [c1, c2] = find(adj);
    
    % Eliminating relations between neighbors that belong to different
    % regions, if requested
    if exist('regions')
        [c1, c2] = restrict_to_regions(c1, c2, regions);
    end
    
    % Computing slope measure between neighbors (as measure, we use the
    % signed, squared tangent of the angle with the horizontal plane)
    neigh_dz  = xyz(c2,3) - xyz(c1,3);
    neigh_dz  = sign(neigh_dz) .* (neigh_dz.^2); % square value, but keep sign
    neigh_dxy = xyz(c2,:) - xyz(c1,:);
    neigh_dxy = sum(neigh_dxy.^2, 2); % the squared xy-distance bwt. neighs
    slopes    = neigh_dz ./ neigh_dxy; % squared tan of angle (or NaN, for
                                       % self-neighbors) 
    slopes(find(isnan(slopes))) = 0; % setting NaN values to zero.
    
    % Assembling the neighborhood matrix
    [dummy, n] = rlencode(c2);
    npos       = cumsum([1;n]);
    nhood      = expand_neighbor_array(c1, npos, 0);

    % Identifying the max. upslope neighbor of each cell (or node)
    slopemat = expand_neighbor_array(slopes, npos, NaN);
    [v, i]   = max(slopemat, [], 2);
    
    upneigh = nhood(sub2ind(size(nhood), 1:num, i'))';
    
end

% ----------------------------------------------------------------------------
function [c1, c2] = restrict_to_regions(c1, c2, regions)

    r1 = regions(c1);
    r2 = regions(c2);
    keep = find(r1 == r2);
    c1 = c1(keep);
    c2 = c2(keep);
end

% ----------------------------------------------------------------------------
function res = expand_neighbor_array(n_ix, npos, initval)
    
    num_cells     = numel(npos) - 1;
    num_neighs    = diff(npos);
    max_num_neigh = max(diff(npos));
    empties       = max_num_neigh - num_neighs;
    intervals     = [num_neighs, empties]';
    vals          = repmat([true;false], num_cells,1);
    indices       = rldecode(vals, intervals);    
    res           = ones(numel(indices),1) * initval;
    res(indices)  = n_ix;
    res           = reshape(res, max_num_neigh, num_cells)';
end
