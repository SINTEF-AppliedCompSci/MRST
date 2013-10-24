function [down_cells,rest_volume,down_cells_prev, Gtop]=spillPointEstimate(Gtop,cell,volume,varargin)
% Find the path of a given volume by spill point analyise
%
% SYNOPSIS:
%     [down_cells,rest_volume,down_cells_prev, Gtop]=spillPointEstimate(Gtop,cell,volume,varargin)
%
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%            supported options are:
%
%              use_multipoint -- locical to extend neighbours node neigbours 
%
% PARAMETERS:
%   G     - Valid top surface grid
%
% RETURNS:
% trap_stuct  - structure describing the traps of the top surface by
%
%   down_cells - cells which can be reached
%   rest_volume - volume left typically negative value
%   down_cells_prev - cells reached by fully filling the level below
%   Gtop - the flat surface corresponding to the spillpoint which will be
%           past with more volume injected
%
%
%  NOTE:
%    
%
%  Example:
%    
%    
%
% SEE ALSO:
%    findTrappingStructure
% 
opt = struct('use_multipoint',false,'include_bc',false);
opt = merge_options(opt, varargin{:});
%require gridtools mex/matlab_bgl
Gtop_new=Gtop;
Gtop_prev=Gtop;
[C,N,CC]=maxTPFAGravityMatrix(Gtop,'use_multipoint',opt.use_multipoint);%#ok
rest_volume=volume;
down_cells_prev=[];
while rest_volume>0,
    % find component of cell
    comp_cells=bfs(C'+C,cell);
    
    
    comp_cells=find(comp_cells>-1);
    % find max of boundary of component         
    if(opt.use_multipoint)
        bc = boundaryCellsSubGrid(Gtop, m,'neigbours',N); 
    else
        bc = boundaryCellsSubGrid(Gtop, comp_cells);
    end
    % flatten top
    [z_level,j]=min(Gtop_new.cells.z(bc));
    % include spill cell which may be in other domain into region
    %%{
    faces=Gtop.cells.faces(Gtop.cells.facePos(bc(j)):Gtop.cells.facePos(bc(j)+1)-1,1);
    faces=unique(faces);
    sp_cells=Gtop.faces.neighbors(faces,:);
    sp_cells=reshape(sp_cells,[],1);
    sp_cells=sp_cells(sp_cells>0);
    sp_cells=unique(sp_cells(:));
    sp_cells=sp_cells(Gtop_new.cells.z(sp_cells)>z_level);
    [z_level_new,j]=min(Gtop_new.cells.z(sp_cells));
    
    if(~isempty(j))
       assert(z_level_new>z_level)
       z_level=z_level_new;
    end
    %}
    figure(33),clf,
    plotGrid(Gtop_new,'FaceColor','none')
    plotGrid(Gtop_new,comp_cells)
    plotGrid(Gtop_new,comp_cells(Gtop_prev.cells.z(comp_cells)<=z_level),'FaceColor','r')
    %plotGrid(Gtop_new,sp_cells,'FaceColor','b')

    comp_cells=comp_cells(Gtop_prev.cells.z(comp_cells)<=z_level);
    Gtop_new.cells.z(comp_cells)=z_level;
    
    % fine new matrix connection
    [C,N,CC]=maxTPFAGravityMatrix(Gtop_new,'use_multipoint',opt.use_multipoint);
    % find cell down stream
    down_cells=bfs(CC,cell);
    if(down_cells(end)>-1)
        % if at boundary
        down_cells=find(down_cells(1:end-1)>-1);
        break;
    end
    
    % calculated stored volume
    down_cells=find(down_cells>-1);
    
    %%{
    z_spill_loc=Gtop_new.cells.z;
z_spill_loc(Gtop.cells.z==Gtop_new.cells.z)=0;
cells=find(z_spill_loc>0);
% find all nodes for this cells
eIX = Gtop.cells.facePos;
nn  = double(diff([Gtop.cells.facePos(cells), ...
   Gtop.cells.facePos(cells + 1)], [], 2));
if(~isfield(Gtop.cells,'cellNodes'))
 Gtop.cells.cellNodes = getCellNodes(Gtop);
end
cn  = double(Gtop.cells.cellNodes(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));
% set all the node values for this cells to the cell center value
zz=rldecode(z_spill_loc(z_spill_loc>0),nn);
Gtop_new.nodes.z(cn)=zz;
%}    
    plotGrid(Gtop_new,'FaceColor','none')
    plotGrid(Gtop_new,down_cells)
    plotGrid(Gtop_new,comp_cells,'FaceColor','r')
    %rest_volume    

    vv=(Gtop_new.cells.z(down_cells)-Gtop_prev.cells.z(down_cells))...
        .*Gtop.cells.volumes(down_cells);
    diff_volume=sum(vv);
    rest_volume=rest_volume-diff_volume;
    Gtop_prev=Gtop_new;
    down_cells_prev=down_cells;
end
z_spill_loc=Gtop_new.cells.z;
z_spill_loc(Gtop.cells.z==Gtop_new.cells.z)=0;
cells=find(z_spill_loc>0);
% find all nodes for this cells
eIX = Gtop.cells.facePos;
nn  = double(diff([Gtop.cells.facePos(cells), ...
   Gtop.cells.facePos(cells + 1)], [], 2));
if(~isfield(Gtop.cells,'cellNodes'))
 Gtop.cells.cellNodes = getCellNodes(Gtop);
end
cn  = double(Gtop.cells.cellNodes(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));
% set all the node values for this cells to the cell center value
zz=rldecode(z_spill_loc(z_spill_loc>0),nn);
Gtop_new.nodes.z(cn)=zz;

figure(33),clf,
plotGrid(Gtop_new,'FaceColor','none')
plotGrid(Gtop_new,down_cells)
plotGrid(Gtop_new,comp_cells,'FaceColor','r')
