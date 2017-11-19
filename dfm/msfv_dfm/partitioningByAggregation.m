function p = partitioningByAggregation(g,dual)
% Create primal partitioning by aggregation with node cells as seeds.
%
% SYNOPSIS:
% p = partitioningByAggregation(g_fine,dual)
%
% PARAMETERS:
%       g - Mrst grid structure
%       dual - dual grid as used in the MSFVM module
%
% OUTPUT:
%       p - The partition vector
%
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.

% create a cell2cell map
isintern = all(g.faces.neighbors~=0,2);

n1 = double(g.faces.neighbors(isintern,1));
n2 = double(g.faces.neighbors(isintern,2));
ind = find(isintern);

if isfield(g.cells,'neighbors')
    n1 = [n1 ; g.cells.neighbors(:,1)];
    n2 = [n2 ; g.cells.neighbors(:,2)];
    ind = [ind ; ind(end)+ (1:size(g.cells.neighbors,1))'];
end
cell2cell = sparse([n1 n2],[n2 n1],[ind ind],g.cells.num,g.cells.num);

% Store the dual grid in a cell
d_cells = cell(3,1);
d_cells{1} = dual.nn;
d_cells{2} = dual.ee;
d_cells{3} = dual.ii;

% get the partitioning
p = getPrimalCells(cell2cell,d_cells);

end

function p = getPrimalCells(cell2cell,d_cells)
% get the partitioning


% the seeds are the node cells
primal_cells_this = d_cells{1};

% the new cells index
primal_ind_this = (1:numel(primal_cells_this))';

% final storage
primal_cells = primal_cells_this;
primal_ind = primal_ind_this;

% we do the aggregation level wise starting with the edge cells
for iter = 2:3

    % get the cells to look within
    if iter == 3 && ~isempty(cells)
         cells = [cells;d_cells{iter}];
    else
         cells = d_cells{iter};
    end

    % aggregate until there are no more cells at this level
    while ~isempty(primal_cells_this)
        % one step correspondence to finding the neighboring cells among the
        % remaing cells
        [primal_cells_this,primal_ind_this,cells] = getTheres(cell2cell, primal_cells_this,primal_ind_this,cells);

        % update the storage
        primal_cells = [primal_cells;primal_cells_this];
        if ~isempty(primal_ind_this)
            primal_ind = [primal_ind;primal_ind_this(:)];
        end

    end
    % make ready for next level
    primal_ind_this = primal_ind;
    primal_cells_this = primal_cells;

end

% sort the indices and store the partitioning.
[~,map] = sort(primal_cells);
p = primal_ind(map);

end

function [there,ind,cells] = getTheres(c2c,here,ind,cells)
% get the neighboring cells among the cells

% find the cell neighbors of th here cells among the remaining cells
[ind_tmp,there_tmp] = find(c2c(here,cells));
if isempty(ind_tmp)
    there = [];
    ind = [];
    return
end

% update the indices
numCells = accumarray(ind_tmp,1);
[~,I] = sort(numCells(ind_tmp));
ind_tmp = ind_tmp(I);
there_tmp = there_tmp(I);

% make sure that every cells is counted only once.
% unique pick out the first occurrence of the cells,
% a fairer distribution may be better
[there_tmp,I] = unique(there_tmp,'first');
ind_tmp = ind_tmp(I);

% map to original cells and indexes
there = cells(there_tmp);
ind = ind(ind_tmp);

% remove the added cells from the cells array
cells = setdiff(cells,there);
end

