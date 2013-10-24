function [dstr_neigh, region, spill_edges] = nodeSpillField(Gt)
% Compute the spill field for the 2D grid 'Gt'.
% The spill field consists of a classification of each of the nodes
% in the grid as belonging to a specific spill region.  A spill region is defined by 
% the nodes that all leads (i.e. guides the flow) up to a local maximum of the z-field. 
% The node at the maximum is flagged as a "sommet" node.
%
% SYNOPSIS
%   [downstream_neigh, region, spill_edges] = spill_field(Gt);
%
% PARAMETERS:
%      Gt - a 2D grid with an associated z-field
%
% RETURNS:    
%   dstr_neigh  - One element per node in Gt, containing the index of the
%                 'most (steepest) downstream' of its neighbor nodes, or '0'
%                 if it is a boundary node or sommet node (no downstream neighbor)
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
    % Establishing neighbor relations between nodes (building sparse incidence matrix)
    node_neighbors = node_incidence_matrix(Gt); % one row and one column per node
    num_nodes = size(node_neighbors, 1);
    
    % Determining which nodes lie on the boundary (these are automatically considered as 
    % belonging to the 'exterior' spill region.
    boundary_node = identify_boundary_nodes(node_neighbors); % one entry per node

    % Initializing vector storing information about regions and downstream neighbors.
    % The real values will be recursively determined later by calling the nested
    % function 'establish_region'.  (NB: 'dstr_neigh', 'region' and region_count are the
    % ONLY variables that will be modified by the nested functions below).
    
    dstr_neigh = repmat(-1, num_nodes, 1);  % one entry per node (-1 means 'uninitialized')
    region = repmat(-1, num_nodes, 1);     % Possible flags will be:
                                           % -1 : uninitializied
                                           % -2 : pending (temporary value used by recursive 
                                           %      algorithm)
                                           %  0 : part of the 'exit' region
                                           % >1 : belongs to an interior region
    region(boundary_node) = 0;             % All boundary nodes considered to be part of the
                                           % 'exit' region. 
    dstr_neigh(boundary_node) = 0;         % We do not track downstream  neighbors for boundary
                                           % nodes.
    region_count = 0;                      % Counting interior regions as they are discovered
    
    % Identifying regions and sommets
    for i = 1:num_nodes
        if region(i) == -1
            % region/downstream neighbor info for this node has not yet been set.  Let's do it here.
            assign_region(i);
        end
    end
    
    % Computing spill edges (see documentation of 'spill_edges' at the top of this function).
    spill_edges = find_spill_edges(Gt, region);
    
    %% --- Recursive helper function requiring access to 'node_info', therefore nested 
    %% --- This function assigns a region number to the node with index 'node_ix'.
    %% --- It might have to recursively determine the region number of neighboring indices
    %% --- before deciding the region number of 'node_ix' itself.
    %% --- This function might modify the following variables in its calling context:
    %% --- * 'dstr_neigh'
    %% --- * 'region'
    %% --- * 'region_count'
    function assign_region(node_ix)

    % Determining region and whether node is a sommet or has an usptream neighbor
        if region(node_ix) >= 0
            % region has already been initialized - nothing to do!
            return
        end

        % Finding the node in the neighborhood, and checking if it is strictly greater
        % than the current node
        [dstr_neigh(node_ix), strictly_higher] = highest_in_neighborhood(node_ix);
        if dstr_neigh(node_ix) == node_ix 
            % The current node is strictly the highest in its neighborhood
            dstr_neigh(node_ix) = 0; % We flag this exlicitly as a sommet node.
            region_count    = region_count+1;
            region(node_ix) = region_count;
        elseif strictly_higher
            % This node belongs to the same region as its strictly higher
            % neighbor
            assign_region(dstr_neigh(node_ix));
            region(node_ix) = region(dstr_neigh(node_ix));  
        else 
            % this node has a neighbour node at exactly the same height.  Keep looking
            % for a higher one.
            region(node_ix) = -2; % we flag our current node as 'pending' with '-1'.
            assign_region(dstr_neigh(node_ix));
            region(node_ix) = region(dstr_neigh(node_ix)); % Node will belong to same region as its
                                                           % neighbor with the same height.
        end
    end % --- End of recursive helper function 'assign_region' ---

    
    %% --- Helper function to 'assign_region' above.  Determines which node is 'highest' (in terms of 
    %% --- z-value) among the five nodes consisting of 'n_ix' and its (up to) four neighbor nodes.
    %% --- This function does not modify any of the variables in its calling context.
    function [hn, strictly_higher] = highest_in_neighborhood(n_ix)

        % looking up all neighbors of this node
        neighs = find(node_neighbors(:, n_ix));

        % ignoring (eliminating) neighbors flagged as 'pending'
        neighs = neighs(find(region(neighs) ~= -2));
        
        % Finding highest neighbor (i.e. the one located at the lowest depth)
        dist = sqrt(sum(bsxfun(@minus,Gt.nodes.coords(neighs,:),Gt.nodes.coords(n_ix,:)).^2,2));
        % to avoid dividing by zero
        dist(dist==0)=max(dist);
        neighs_depth = Gt.nodes.z(neighs);
        cur_depth = Gt.nodes.z(n_ix);        
        %[min_neigh_depth, min_neigh_ix] = min(neighs_depth);
        [min_neigh_depth, min_neigh_ix] = min(bsxfun(@minus,neighs_depth,cur_depth)./dist);%#ok
        min_neigh_depth=neighs_depth(min_neigh_ix);
        % What is the depth of the current node?
        

        % Filling in return variables (index of highest node in neighborhood and
        % whether it is strictly higher than the current node)
        strictly_higher = (min_neigh_depth < cur_depth);
        if isempty(min_neigh_depth) || min_neigh_depth > cur_depth
            hn = n_ix;
        else
            hn = neighs(min_neigh_ix);
        end
    end
end

function imat = node_incidence_matrix(Gt)
% Establish the matrix of incidences between nodes in Gt.  The (sparse) matrix has one
% line and one column per node.  Nonzero elements are those (i,j) where node i an node
% j share an edge.

    % Since we are working on a 2D grid, we assume that all edges have two distinct
    % nodes, so we do not bother using the 'Gt.faces.nodePos' indirection map.

    node_incidences = reshape(Gt.faces.nodes, 2, [])';

    nnum = Gt.nodes.num;    
 
    % mapping incidences.  Since, if (i,j) is an incidence, (j,i) must be one too, we
    % make the matrix symmetric in the second step.
    % imat = accumarray(node_incidences, 1, [nnum, nnum], @prod, [], true); % doesn't work in Octave
    tmp = accumarray(node_incidences, ones(size(node_incidences, 1), 1), [nnum, nnum], @prod, [], true);
    tmp = tmp + tmp';

    % EXPLANATION OF THE UPCOMING CODE LINE:
    % Currently, 'tmp' refers to the four nodes reachable from the node itself.  We want to also 
    % include the four (for interior nodes) diagonally reachable nodes in the neighborhood.  These
    % can be identified as those reachable in two steps in at least two ways from the original node,
    % (moving along connecting edges).  We compute this below, and 'OR' it
    % with the original incidence matrix.
    % tmp*tmp gives all nodes reachable in two steps, whereas the test '>1' restricts the result to
    % those reachable in at least two different ways.  (This has also the side effect of including
    % the original node in the neighborhood, but this shouldn't be a problem here).
    imat = tmp | (tmp*tmp>1);

end

function bnodes = identify_boundary_nodes(incidence_matrix)

% Return a vector with one entry per global node.  The entry is 'true' for boundary
% nodes, and 'false' for interior nodes.  A node is considered on the boundary if it
% has less than 9 neighbors (as interior nodes are expected to have, counting itself).  
% There might be some cases where this does not hold exactly, so a more rigorous 
% (but likely slower) implementation might be considered in the future (based on 
% explicit topology information from the grid itself).
    
    bnodes = sum(incidence_matrix) < 9;
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
    % consecutively in the Gt.faces.nodes vector (@@ is this always the case?)
    enodes = reshape(Gt.faces.nodes, 2, [])';
    
    %      --- Col 1 ---     -- Col 2 and 3 -- -- Col 4 and 5 --  -- Col 6 and 7 --
    res = [(1:Gt.faces.num)', s_field(enodes), Gt.nodes.z(enodes), double(enodes)];
    % NB: 'enodes' converted to 'double' in column 6 and 7 above to avoid
    % implicit conversion of the whole 'res' matrix into 'integer' (this
    % initially created a subtle, difficult-to-track-down bug).
    
    % Keep only the lines where the edge's two nodes belong to different spill regions
    res = res(res(:,2) ~= res(:,3), :);

end












