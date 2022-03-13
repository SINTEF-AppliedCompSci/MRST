function [uslp_neigh, region, spill_edges] = nodeSpillField(Gt, closed_bnodes, closed_fnodes)
% Compute the spill field for the 2D grid 'Gt'.
% The spill field consists of a classification of each of the nodes
% in the grid as belonging to a specific spill region.  A spill region is defined by 
% the nodes that all leads (i.e. guides the flow) up to a local maximum of the z-field. 
% The node at the maximum is flagged as a "sommet" node.
%
% SYNOPSIS
%   [upslope_neigh, region, spill_edges] = spill_field(Gt);
%
% PARAMETERS:
%    Gt            - a 2D grid with an associated z-field
%    closed_bnodes - vector with indices of nodes on the closed part of the
%                    boundary 
%    closed_fnodes - vector with indices of nodes along closed fault lines
% 
% RETURNS:    
%   uslp_neigh  - One element per node in Gt, containing the index of the
%                 'most (steepest) upslope' of its neighbor nodes, or '0'
%                 if it is a boundary node or sommet node (no upslope neighbor)
%   region      - a vector with one element per node in Gt.  The element gives the index
%                 of the spill region that the corresponding node belongs to.  Nodes with
%                 index 0 belong to the 'exterior' spill region.  Nodes with higher indices
%                 belong to a specific interior spill region.
%   spill_edges - Information about all spill edges in the grid.  A spill
%                 edge is defined as an edge whos two nodes belong to
%                 different spill region.  'spill_edges' is a matrix with one
%                 line per identified spill edge and seven columns, containing
%                 the following information about the spill edge:
%                    Col 1: index of spill edge (into Gt.faces)
%                    Col 2: region of edge's first node
%                    Col 3: region of edge's second node
%                    Col 4: z-value of edge's first node
%                    Col 5: z-value of edge's second node
%                    Col 6: index of edge's first node
%                    Col 7: index of edge's second node
%
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

    %%% Determining the neighborhood and steepest upslope neighbor of each node
    % Establishing list of all neighborhood relations between nodes.  Since
    % we are working on a 2D grid, we assume that all edges have two distinct
    % nodes, so we do not bother to use the 'Gt.faces.nodePos' indirection map.
    neigh_rels      = reshape(Gt.faces.nodes, 2, [])';
    fnode_indicator = sum(ismember(neigh_rels, closed_fnodes(:)) > 0, 2);
    neigh_rels      = neigh_rels(~fnode_indicator,:); % remove relations
                                                      % involving nodes of
                                                      % closed faults
    nodes_xyz       = [Gt.nodes.coords, Gt.nodes.z];
    [uslp_neigh, nhood] = findUpslopeNeighbor(nodes_xyz, neigh_rels);
                                               
    
    %%% Constructing the 'region' matrix
    % Flag boundary nodes as those having less than nine members of
    % their neighborhood, causing (at least) the last index of 'nhood' to be
    % zero. 
    %boundary_nodes = (nhood(:,end) == 0);
    
    boundary_edges = find(prod(Gt.faces.neighbors,2)==0);
    boundary_nodes = false(Gt.nodes.num, 1);
    boundary_nodes(Gt.faces.nodes(mcolon(Gt.faces.nodePos(boundary_edges), ...
                                         Gt.faces.nodePos(boundary_edges+1)-1))) = true;
    open_boundary_nodes = boundary_nodes;
    open_boundary_nodes(closed_bnodes) = false;

    int_sommets            = (uslp_neigh == (1:Gt.nodes.num)') & ~open_boundary_nodes;
    region                 = NaN * ones(Gt.nodes.num,1);
    region(open_boundary_nodes) = 0;
    region(int_sommets)    = 1:sum(int_sommets); % assigning region identifiers
    unfinished             = find(isnan(region));

    if ~all(boundary_nodes(closed_bnodes))
       warning(['Ignoring some nodes specified as closed boundary nodes, as ' ...
                'they were not found on the boundary.']);
    end

    while ~isempty(unfinished)
        cur_sneigh_ix = uslp_neigh(unfinished);
        defineds = find(~isnan(region(cur_sneigh_ix)));
        actives = unfinished(defineds);
        region(actives) = region(cur_sneigh_ix(defineds));
        
        unfinished = find(isnan(region));    
    end
    
    %%% Final adjustments of uslp_neigh, to fulfill contract
    % Ascribe boundary nodes and interior sommet nodes an upslope neighbor
    % index of 0
    uslp_neigh(open_boundary_nodes) = 0;
    uslp_neigh(int_sommets) = 0;

    % Computing spill edges
    spill_edges = find_spill_edges(Gt, region);
end

function res = find_spill_edges(Gt, s_field)
% Determine which edges in the 2D grid are 'spill_edges'.  A spill edge is here defined
% as an edge whose two endpoint nodes belong to separate spill regions.  
% 
% SYNOPSIS
%   res = find_spill_edges(Gt, spill_field);
%
% PARAMETERS:
%            Gt - a 2D grid with assocated z-field.  The function will identify the
%                 spill edges of this grid.
%   spill_field - The precalculated spill_field of Gt (as calculated within the 
%                 'spill_field' function.
%
% RETURNS:
%   res - a matrix with seven columns and one line per identified spill edge.  The
%         columns represent the following information:
%         Col 1: index of spill edge (into Gt.faces)
%         Col 2: region of edge's first node
%         Col 3: region of edge's second node
%         Col 4: z-value of edge's first node
%         Col 5: z-value of edge's second node
%         Col 6: index of edge's first node
%         Col 7: index of edge's second node
    
    % We here assume that each face (i.e. edge) has exactly two nodes, which are stored 
    % consecutively in the Gt.faces.nodes vector 
    enodes = reshape(Gt.faces.nodes, 2, [])';
    
    %      --- Col 1 ---     -- Col 2 and 3 -- -- Col 4 and 5 --  -- Col 6 and 7 --
    res = [(1:Gt.faces.num)', s_field(enodes), Gt.nodes.z(enodes), double(enodes)];
    % NB: 'enodes' converted to 'double' in column 6 and 7 above to avoid
    % implicit conversion of the whole 'res' matrix into 'integer' (this
    % initially created a subtle, difficult-to-track-down bug).
    
    % Keep only the lines where the edge's two nodes belong to different spill regions
    res = res(res(:,2) ~= res(:,3), :);

end

