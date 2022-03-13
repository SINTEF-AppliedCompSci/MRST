function trap=findTrappingStructure(Gt, varargin)
%Find the trapping structure of a top-surface grid
%
% SYNOPSIS:
%   trap=findTrappingStructure(Gt)
%
% PARAMETERS:
%   Gt   - Valid top-surface grid
%
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%            supported options are:
%
%              use_multipoint -- locical to extend neighbours to node based
%              stencil for neighbours.
%
% RETURNS:
%   a structure describing the traps of the top surface. The structure
%   consists of the following fields (nc=number of cells in Gt):
%
%   top  - an (nx1) vector of cell numbers for all local maxima on the top
%          surface grid
%
%   istrap - an (nx1) boolean vector identifying whether a local maximum is
%          a trap or not. In the following, m = numel(top(istrap))
%
%   Gtop - a top-surface grid, in which the topography inside each
%          (primary) trap has been replaced by a flat surface at the depth
%          of the corresponding spill point
%
%   z_spill_loc - an (nc x 1) vector giving the spill-point depth of the
%          largest trap a given cell is part of, and zero if the cell is
%          not part of a trap.
%
%   g_trap - a (nc x m) sparse matrix. Each column corresponds to a local
%          maximum and has nonzero entries for all cells that are part of
%          the same trap as the local maximum point. Positive values give
%          the maximum 'trapping level'. A level-1 trap contains one local
%          maximum point, a level-2 trap contains at least two level-1
%          traps and hence at least one local spill point, a level-3 trap
%          contains at least two level-2 traps, and so on.
%
%   trap_level - a (tl x 1) cell array of (nc x m) sparse matrices, where
%          'tl' is the maximum trap level. Each sparse matrix describes the
%          traps on a given level, that is, trap_level{k}{i,j} is equal one
%          if cell number i is part of a level-k trap that has cell j as a
%          local maximum.
%
%   z_spill_level - a (tl x 1) cell array of (m x 1) vectors defined so
%          that z_spill_level{k}(j) gives the spill-point depth of the
%          level-k trap that the local maximum in cell j is part of and NaN
%          if cell j is not a local maximum in a level-k trap.
%
%   z_spill_loc_level - a (tl x 1) cell array of (nc x 1) vectors. The
%          vector z_spill_loc_level{k}(i) gives the depth of the spill
%          point for the level-k trap that cell i is part of and zero if
%          cell is not part of a trap.
%
% NOTE:
%   The routine assumes that all top-grid surfaces have strictly positive
%   z-values and will fail if the z-values cross zero.
%
% Examples:
%   ts=findTrappingStructure(Gt)
%   
%   % Plot the new top-surface grid on top of the old one to show
%   % the location of traps
%   plotGrid(Gt,      'FaceColor','r','EdgeAlpha',.1)
%   plotGrid(ts.Gtop, 'FaceColor','y','EdgeAlpha',.1);
%   view(3), axis tight off
%
%   % Show the flat parts of the grid
%   cla
%   plotGrid(ts.Gtop,'FaceColor','none','edgeAlpha',0.1)
%   plotGrid(ts.Gtop, ts.z_spill_loc>0)
%   view(3), axis tight off
%
%   % Show the traps on different levels in different colors
%   cla
%   plotGrid(ts.Gtop,'FaceColor','none','edgeAlpha',0.1)
%   nt  = numel(ts.trap_level);
%   col = hsv(nt);
%   for k=nt:-1:1
%      cell_list = [];
%      for i=1:size(ts.trap_level{k},2)
%         cell_list = [cell_list; find(ts.trap_level{k}(:,i)>0)];
%      end
%      plotGrid(Gt, cell_list, 'FaceColor', col(k,:), 'EdgeAlpha',.1);
%   end
%   colormap(col); cbh=colorbar; caxis([0 nt]+.5);set(cbh,'YTick',1:nt);
%
% SEE ALSO:
%   `topSurfaceGrid`, `trapAnalysis`, `findCellLines`, `findTrapConnections`

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

opt = struct('use_multipoint',false);
opt = merge_options(opt, varargin{:});

mlist = mrstModule();
mrstModule add matlab_bgl;

% Find top cells of the surface defined such that all centroids of all
% internal neighbors lie deeper.
nc       = Gt.cells.num;

if(opt.use_multipoint)   
  N=neighboursByNodes(Gt,'include_bc',true);
else
   N=Gt.faces.neighbors;
end
%
internal=all(N>0,2);
ni1=N(internal,1);ni2=N(internal,2);
z_diff   = Gt.cells.z(ni2)-Gt.cells.z(ni1);
top      = accumarray([ni2;ni1], double([z_diff;-z_diff]<=0), [nc,1], @prod, 0);
%top      = accumarray(ni2, double(z_diff<0), [nc,1], @prod, 0);
top_cells=find(top);
%
% Make copy of input grid which is going to be modifed
Gtop=Gt;
 
% Loop through all cells that contain a local maximum and find the
% corresponding trap through a breadth-first search. If a trap of level>1
% contains another local maximum, this local maximum is removed from the
% list of start positions for subsequent trap searches.
candidate = false(nc,1); candidate(top_cells)=true;
has_trap = false(numel(top_cells),1);
trap_level={};
kk=0;
while any(candidate)
   kk=kk+1;
   z_diff = Gtop.cells.z(ni2) - Gtop.cells.z(ni1);
   A = sparse(double(ni2), double(ni1), z_diff>=0, nc, nc);
   A = A + sparse(double(ni1), double(ni2), z_diff<=0, nc, nc);

   z_spill_level{kk}  = nan(numel(top_cells),1);%#ok
   
   trap_ind=[];   
   for i=1:numel(top_cells)
      if ~candidate(top_cells(i)), continue; end
      
      % Do a breadth-first search to find the trap defined as all cells
      % that are reachable from the local maximum by moving out of a cell
      % and into a neighbouring cells whos centroid lies deeper than the
      % previous cell. The search picks the neighbour with the largest
      % depth difference and terminates when encountering a cell with a
      % centroid that is higher than in the current cell.
      [istrap, z_spill] = bfs_argumented(A,N, Gtop, top_cells(i), opt);      
      
      % extract the cell number of all cells in the trap. If one of the
      % cells contains a local maximum, this cell is removed from the list
      % of start points for subsequent trap searches. Replace the
      % topography inside the trap with a flat surface at the elevation of
      % the spill point
      trap_cells = find(istrap);
      if numel(trap_cells)
         has_trap(i) = true;
         candidate(trap_cells) = false;
         candidate(top_cells(i)) = true;
         trap_ind=[trap_ind;trap_cells,repmat(i,numel(trap_cells),1)];%#ok
 
         Gtop.cells.z(trap_cells)=z_spill;
      else
         candidate(top_cells(i)) = false;
      end
      
      z_spill_level{kk}(i)=z_spill;%#ok
   end
   
   % define sparse matrix that defines all level-kk traps
   if ~isempty(trap_ind)
      fprintf(1,'Trap level %d: %d traps identified\n', kk, sum(candidate));
      trap_level{kk}=sparse(trap_ind(:,1),trap_ind(:,2),1,nc,numel(top_cells));%#ok
   else
      z_spill_level = z_spill_level(1:kk-1);
   end
   
end

% Remove local maxima that have no other cells attached (these typically
% appear at the boundary of the domain).
for i=1:numel(trap_level)
   z_spill_level{i} = z_spill_level{i}(has_trap);%#ok
   trap_level{i} = trap_level{i}(:,has_trap);    %#ok
end

% Check for consistency between different trap levels. Determine the
% spill point depth for all cells inside a trap and the maximum spill-point
% depth for each cell. If a cell is not part of any trap, the value is
% zero. For this to work, the top-surface grid must have positive z-values
% only (i.e., the surface must not cross the plane z=0).
z_spill_loc=zeros(nc,1);
for kk=1:numel(trap_level)
   z_spill_loc_level{kk}=z_spill_loc;%#ok
   for i=1:numel(top_cells(has_trap))
      ind=find(trap_level{kk}(:,i)>0);
      if(~isempty(ind))
         z_spill_loc_level{kk}(ind)=max(z_spill_loc(ind),z_spill_level{kk}(i));%#ok
         z_spill_loc(ind)=max(z_spill_loc(ind),z_spill_level{kk}(i));
      end
   end         
end

% Define a vector with the maximum trap level of each cell 
if(numel(trap_level)>1)
    g_trap=trap_level{1};   
else
    g_trap=[];
    z_spill_level={};
    z_spill_loc_level={};
end
for kk=2:numel(trap_level)
   for i=1:size(trap_level{kk},2)
      % all higher traplevels are included in smaller traps
      ok=trap_level{kk}(trap_level{kk-1}(:,i)>0,i)==1;
      if(~all(ok))
         assert(all(trap_level{kk}(:,i)==0));        
      end
   end
   g_trap(g_trap==0 &trap_level{kk})=kk;
end

% Modify the z-values of the nodes for all cells that are inside a trap
% equal to the z-value of the cell center
cells     = find(z_spill_loc>0); 
eIX       = Gtop.cells.facePos;      
nn        = double(diff([Gtop.cells.facePos(cells), ...
                         Gtop.cells.facePos(cells + 1)], [], 2));


Gtop.cells.sortedCellNodes = getSortedCellNodes(Gtop);

cn        = double(Gtop.cells.sortedCellNodes(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));
zz        = rldecode(z_spill_loc(z_spill_loc>0),nn);
Gtop.nodes.z(cn) = zz;

z_spill_loc(z_spill_loc==Gt.cells.z)=0;
top_cells=top_cells(z_spill_loc(top_cells)>0);

% Pack the results in a data structure
trap =struct('Gtop', Gtop, ...
             'top',  top_cells, ...
             'istrap', has_trap, ...
             'z_spill_loc',z_spill_loc,...
             'g_trap',g_trap,...
             'trap_level',{trap_level},...
             'z_spill_level',{z_spill_level},...
             'z_spill_loc_level',{z_spill_loc_level});
                     
mrstModule('reset', mlist{:})
return
end


function [mycells_p, z_spill]= bfs_argumented(A, N, Gtop, start_cell, opt)
% Breadth first search, designed to first go in the direction of maximal
% difference in height values of the centroids of two neighbouring cells
require('matlab_bgl')
% search all cells backwards connected
cc=dfs(A', start_cell);
cc=cc+1;cc(cc>0)=2;
m=find(cc>0);

% find the cells on the edge of the domain
% equivalent to for normal neighbours
if(opt.use_multipoint)
       bc = boundaryCellsSubGrid(Gtop, m,'neigbours',N); 
else
    bc = boundaryCellsSubGrid(Gtop, m);
end
[z_spill,j]=min(Gtop.cells.z(bc));%#ok
% all cells below or at the same level as the spill point
mycells_p=(cc>0 & Gtop.cells.z<=z_spill);

% all cells above the spill point connected to the start cell,
% define trap as empty if no cells have smaller z-value than
% the spill point
if(~any((cc>0 & Gtop.cells.z<z_spill)))
    mycells_p=zeros(size(mycells_p));
end
end
