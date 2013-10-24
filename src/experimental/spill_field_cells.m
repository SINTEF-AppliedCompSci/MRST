function [sommet, region] = spill_field_cells(Gt)
% Compute the spill field for 2D grid 'Gt'.
% The spill field consists of a classification of each of the cells
% in the grid as belonging to a specific region.  A region is defined by the ensemble
% of cells that all leads (i.e. guides the flow) up to a local maximum of the z-field. 
% The cell at the maximum is flagged as a "sommet" node.
%
% SYNOPSIS
%   [sommet, region] = spill_field(Gt);
%
% PARAMETERS:
%      Gt - a 2D grid with an associated z-field
%
% RETURNS:    
%   sommet - a vector with one element per cell in Gt.  The element is 'true' for all
%            cells identified as sommets (local maxima), and 'false' for all others.
%   region - a vector with one element per cell in Gt.  The element gives the index
%            of the spill region that the corresponding cell belongs to.  Cells with
%            index 1 belong to the 'exterior' spill region.  Cells with higher indices
%            belong to a specific interior spill region.
    
    % Establishing neighbor relations between cells
    internal=all(Gt.faces.neighbors>0,2);      
    cell_neighbors = cell_incidence_matrix(Gt,internal); % one row and one column per cell


    % Determining which cells lie on the boundary    
    cellno=rldecode([1:Gt.cells.num]',diff(Gt.cells.facePos));
    boundary_cells = ~accumarray(cellno,double(internal(Gt.cells.faces(:,1))),[Gt.cells.num,1],@prod,[],true);
    
    % Initializing vector storing information about regions and flagging sommets.
    % The real values will be recursively determined later by calling the nested
    % function 'establish_region'.  (NB: 'sommet', 'region' and region_count are the
    % ONLY variables that will be modified by the nested functions below).
    
    sommet = repmat(false, Gt.cells.num, 1);  % one entry per cell
    region = zeros(Gt.cells.num, 1);          % one entry per cell (0 flags: 'uninitialized')
    region(boundary_cells) = 1;             % Identifying the 'exit' region (which starts on
                                           % the boundary).  The 'exit' region is flagged by
                                           % '1' (interior regions will be flagged with
                                           % numbers from 2 and upwards)
    region_count = 1;                      % We count the exit region as the first region
    
    % Identifying regions and sommets
    for i = 1:Gt.cells.num
        if region(i) == 0
            % region/sommet info for this cell has not yet been set.  Let's do it here.
            assign_region(i);
        end
    end
    
    
    %% Recursive helper function requiring access to cell_info, therefore nested
    function assign_region(cell_ix)

    % Determining region and whether cell is a sommet
        if region(cell_ix) > 0
            % region has already been initialized - nothing to do!
            return
        end

        % Finding the cell in the neighborhood, and checking if it is strictly greater
        % than the current cell
        [hn_ix, strictly_higher] = highest_in_neighborhood(cell_ix);
        assert(hn_ix<Gt.cells.num)
        if hn_ix == cell_ix 
            % The current cell is strictly the highest in its neighborhood
            sommet(cell_ix) = true;
            region_count    = region_count+1;
            region(cell_ix) = region_count;
        elseif strictly_higher
            % This cell belongs to the same region as its strictly higher neighbor
            assign_region(hn_ix);
            region(cell_ix) = region(hn_ix);  % sommet(cell_ix) remains 'false'
        else 
            % this cell has a neighbour cell at exactly the same height.  Keep looking
            % for a higher one.
            region(cell_ix) = -1; % we flag our current cell as 'pending' with '-1'.
            assign_region(hn_ix);
            region(cell_ix) = region(hn_ix); % Cell will belong to same region as its
                                             % neighbor with the same height.
        end
    end 
    %% End of recursive helper function 'assign_region'
    
    %% Helper function to assign_region
    function [hn, strictly_higher] = highest_in_neighborhood(n_ix)

        % looking up all neighbors of this cell
        neighs = find(cell_neighbors(:, n_ix));

        % ignoring (eliminating) neighbors flagged as 'pending'
        neighs = neighs(find(region(neighs) ~= -1));
        
        % Finding highest neighbor (i.e. the one located at the lowest depth)
        neighs_depth = Gt.cells.z(neighs);
        [min_neigh_depth, min_neigh_ix] = min(neighs_depth);
        
        % What is the depth of the current cell?
        cur_depth = Gt.cells.z(n_ix);

        % Filling in return variables (index of highest cell in neighborhood and
        % whether it is strictly higher than the current cell)
        strictly_higher = (min_neigh_depth < cur_depth);
        hn = if_else(isempty(min_neigh_depth) | min_neigh_depth > cur_depth, ...
                     n_ix, ...
                     neighs(min_neigh_ix));
    end
end

function imat = cell_incidence_matrix(Gt,internal)
% Establish the matrix of incidences between nodes in Gt.  The (sparse) matrix has one
% line and one column per node.  Nonzero elements are those (i,j) where node i an node
% j share an edge.

    % Since we are working on a 2D grid, we assume that all edges have two distinct
    % nodes, so we do not bother using the 'Gt.faces.nodePos' indirection map.

    % mapping incidences.  Since, if (i,j) is an incidence, (j,i) must be one too, we
    % make the matrix symmetric in the second step.
    imat=sparse(double(Gt.faces.neighbors(internal,1)),double(Gt.faces.neighbors(internal,2)),1,Gt.cells.num,Gt.cells.num);
    imat=imat+imat';
    
end

function bnodes = identify_boundary_calls(incidence_matrix)

% Return a vector with one entry per global node.  The entry is 'true' for boundary
% nodes, and 'false' for interior nodes.  A node is considered on the boundary if it
% has less than 4 neighbors.  There might be some cases where this does not hold
% exactly, so a more rigorous (but likely slower) implementation might be considered in
% the future (based on additional information found in Gt).
    
    bnodes = sum(incidence_matrix) < 4;
    
end



















