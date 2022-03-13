function trap_con = findTrapConnections(Gnew,z_spill_loc)
% Find the connections between the traps in a top-surface grid
%
% SYNOPSIS:
%   trap_con = findTrapConnections(Gnew,z_spill_loc)
%
% PARAMETERS:
%   Gnew - valid top surface grid 
%
%   z_spill_loc - Vector of size number of cells which define traps. The
%         vector is >0 and equal the z value of the traing surface
%         for cells in a trap and =0 otherwise
%
% RETURNS:
%   trap_con - A struct with connection information between the traps.
%               It has the form
%
%        trap_con=struct(  'cell_lines',cell_lines,...
%                          'cell_line_traps', cell_line_traps, ...
%                          'traps',traps,...
%                          'trap_matrix',trap_matrix,...
%                          'leaf_lines',leaf_lines,...
%                          'leaf_traps',leaf_traps);
%
%       cell_lines - Cell_array of arrays containing the cells connectin
%                    two traps including the spill cell and entry cell
%       cell_line_traps - (l, 2)-matrix, where 'l'=number of entries in 
%                    'cell_lines'.  The two columns of cell_line_traps give
%                    the indices of  the start and destination trap for 
%                    the corresponding cell line.
%
%       traps      - A vector size cells where the number indicate which
%                    trap the cell is assosiated with 0 indicate no trap
%                    number_of_traps=max(traps);
%
%       trap_matrix - Matrix size(number_of_traps,number_of_traps)
%                     describing the acyclic graph connecting the traps
%
%       leaf_lines - Cellarray of vectors with all cells of which connect
%                    the traps upstream of the leaftrap, the line
%                    assosiated with the lowest leafnode first.
%    
%       leaf_traps - Trapnumber of each leafnode, one for each leaf_line and
%                    and sorted by lowest first
%
%
% Example:
%     trap_str = findTappingStructure(Gt) 
%     trap_con = findTrapConnections(Gt,z_spill_loc)
%     figure(),clf,
%     plotGrid(Gnew,'FaceColor','none','edgeAlpha',0.1)
%     plotGrid(Gnew, trap_str.z_spill_loc>0),axis off;axis tight;view(3)  
%     plotGrid(Gnew,{trap_con.cell_lines{:}})
%     axis tight; axis off;view(2);
% 
% SEE ALSO:
%   `topSurfaceGrid`, `trapAnalysis`, `findTrappingStructure`

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

mlist = mrstModule();
mrstModule add coarsegrid matlab_bgl


% make top surface flat according to z_spill_loc
cells=find(z_spill_loc>0);
Gnew.cells.z(cells)=z_spill_loc(cells);
eIX = Gnew.cells.facePos;
nn  = double(diff([Gnew.cells.facePos(cells), ...
   Gnew.cells.facePos(cells + 1)], [], 2));
cellNodes = getSortedCellNodes(Gnew);
cn  = double(cellNodes(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));
zz=rldecode(z_spill_loc(z_spill_loc>0),nn);
Gnew.nodes.z(cn)=zz;
Gnew.cells.cellNodes = getSortedCellNodes(Gnew);

% find connected traps

z_unique = unique(z_spill_loc);
%assert(z_unique(1) == 0);
traps = zeros(size(z_spill_loc));
for i = 1:size(z_unique, 1)
    traps(z_spill_loc == z_unique(i)) = i;
end


% this need the mrst coarse grid module
traps = processPartition(Gnew,traps);
traps = traps-1;

% purge accidentally created traps in non-trap regions 
bad_traps = sort(unique(traps(intersect(find(z_spill_loc == 0), find(traps)))), ...
                 'descend');
if ~isempty(bad_traps)
    fprintf('Warning: purging %d bad traps with 0 spill value.\n', ...
            numel(bad_traps));
    for bt = bad_traps'
        traps(traps == bt) = 0;
        traps(traps > bt) = traps(traps > bt) - 1;
    end
end

I=[];J=[];
traps_z=nan(max(traps),1);
for i=1:max(traps)
   cells=find(traps==i);
   traps_z(i)=Gnew.cells.z(cells(1));
end

% inizialise variables to be filled
vizited_traps=zeros(max(traps),1); % no traps visited
%cell_line=[];
cell_lines={};
leaf_lines={};
leaf_traps={};
count=1;% conting riveres connecting traps
leaf_count=1;% counting leaf nodes
disp('Start find rivers')
while any(vizited_traps==0)
   nv=find(vizited_traps==0);
   % find the lowest trap
   [zz,kk]=max(traps_z(nv));%#ok
   
   trap=nv(kk);
   tstop=false;
   cell_line=[]; % will contain all cellls in a line upstream if a leaf node   
   leaf_traps{leaf_count}=trap;%#ok
   while ~tstop
      % trap visited stop
      if( vizited_traps(trap) )         
         % tstop=true;
         break;%??
      end
      vizited_traps(trap)=true;
      cells=find(traps==trap);
      % simple search in grid (could have used processed partitions)
      % find outflow cell
      cf=Gnew.cells.faces(mcolon(Gnew.cells.facePos(cells),Gnew.cells.facePos(cells+1)-1),1);
      cc=Gnew.faces.neighbors(cf,:);
      % lines below indicate some trouble with the  trapse
      cc=cc(:);
      cc=unique(cc);
      cc=setdiff(cc,[cells;0]);

      uc=find(Gnew.cells.z(cc)<=Gnew.cells.z(cells(1)));
      % if trap neibours boundary stop
      if(numel(uc)==0)
         %tstop=true;
         break;
      else
          % find neigbourng cell with minium z value
         [z,kk]=min(Gnew.cells.z(cc(uc)));%#ok
         uc=uc(kk);
      end
      % cell with first flow
      cell=cc(uc);  
      % cell with first flow 
      cell_line=[cell_line,cell];%#ok
      cell_lines{count}=cell;%#ok
      I=[I;trap];%#ok
      % finding the line connecting the trap staring from cell
      % to the next one
      stopp=false;
      while ~stopp
         % find outflow cells 
         cf=Gnew.cells.faces(Gnew.cells.facePos(cell):Gnew.cells.facePos(cell+1)-1,1);
         cc=Gnew.faces.neighbors(cf,:);
         cc=cc(:);
         cc=unique(cc);
         % if reached boundary stop
         if(any(cc==0))
            stopp=true;
         else
            [z,k]=min(Gnew.cells.z(cc));            
            if(z>=Gnew.cells.z(cell))
               stopp=true;
            else
               % stop if we reach a trap
               cell=cc(k);
               if(traps(cell)>0)
                   stopp=true;                                             %#ok
                   break;
               end
               cell_line=[cell_line,cc(k)];%#ok
               cell_lines{count}=[cell_lines{count},cc(k)];%#ok
               
            end
         end
      end
      count=count+1;
      trap=traps(cell);
      J=[J;trap];%#ok
      % it may be that we are not in a trap cell we should be at boundary
      if(trap==0)
         tstop=true;
      end
   end
   leaf_lines{leaf_count}=cell_line;%#ok
   leaf_count=leaf_count+1;
end
clinetraps = [I, J];
I=I(J>0);
J=J(J>0);
trap_matrix=sparse(I,J,1,max(traps),max(traps));

trap_con=struct('cell_lines',{cell_lines},...
                'cell_line_traps', clinetraps, ...
                'traps',traps,...
                'trap_matrix',trap_matrix,...
                'leaf_lines',{leaf_lines},...
                'leaf_traps',{leaf_traps});
mrstModule('reset', mlist{:})
return;
end
