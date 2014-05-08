function [cell_lines,traps,trap_matrix,leaf_lines,leaf_traps]=findCellLines(Gnew,z_spill_loc)
cells=find(z_spill_loc>0);
eIX = Gnew.cells.facePos;
nn  = double(diff([Gnew.cells.facePos(cells), ...
   Gnew.cells.facePos(cells + 1)], [], 2));
cellNodes = getCellNodes(Gnew);
cn  = double(cellNodes(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));
zz=rldecode(z_spill_loc(z_spill_loc>0),nn);
Gnew.nodes.z(cn)=zz;
Gnew.cells.cellNodes=getCellNodes(Gnew);
%GG=computeGeometryVE(Gnew)
traps=zeros(Gnew.cells.num,1);

traps=double(z_spill_loc>0)+1;
% this need the mrst coarse grid module
require coarsegrid
traps = processPartition(Gnew,traps);
traps = traps-1;
%for i=1:max(p);
%   
%end
   
I=[];J=[];
% number the different traps
% sems to work best fo  some reason
%%{
%for i=1:size(g_trap,2)
%   traps(g_trap(:,i)>0)=i;
%end
%}
%traps_name=unique(traps(traps>0));
for i=1:max(traps)
   cells=find(traps==i);
   traps_z(i)=Gnew.cells.z(cells(1));
end

%%
%trap_mat=sparse(max(traps),max(traps))
%vizited_traps=zeros(numel(traps_name),1);
vizited_traps=zeros(max(traps),1);
%cell_line=[];
cell_lines={};
leaf_lines={};
leaf_traps={};
count=1;
leaf_count=1;
disp('Start find rivers')
tic;
while any(vizited_traps==0)
   nv=find(vizited_traps==0);
   [zz,kk]=max(traps_z(nv));
   %trap=traps_name(nv(kk));
   trap=nv(kk);
   %trap=3;
   tstop=false;
   cell_line=[];
   leaf_traps{leaf_count}=trap;
   while ~tstop       
      %vizited_traps(find(trap==traps_name))=true;
      if( vizited_traps(trap) )
         tstop=true;
      end
      vizited_traps(trap)=true;
      cells=find(traps==trap);
      % simple search in grid
      % find outflow cell
      cf=Gnew.cells.faces(mcolon(Gnew.cells.facePos(cells),Gnew.cells.facePos(cells+1)-1),1);
      cc=Gnew.faces.neighbors(cf,:);
      % lines below indicate some trouble with the  trapse
      cc=cc(:);
      cc=unique(cc);
      cc=setdiff(cc,[cells;0]);
      %[z,kk]=min(Gnew.cells.z(cc));
      %uc=find(Gnew.cells.z(cc)<=Gnew.cells.z(cells(1)));
      uc=find(Gnew.cells.z(cc)<=Gnew.cells.z(cells(1)));
      %assert(numel(uc)>=1);
      if(numel(uc)==0)
         tstop=true;
         break;
      else
         [z,kk]=min(Gnew.cells.z(cc(uc)));
         uc=uc(kk);
      end
      cell=cc(uc);
      %cell=cells(1)
      cell_line=[cell_line,cell];
      cell_lines{count}=cell;
      I=[I;trap];
      % finding the line conecting the trap staring from cell
      % to the next one
      stopp=false;
      while ~stopp
         cf=Gnew.cells.faces(Gnew.cells.facePos(cell):Gnew.cells.facePos(cell+1)-1,1);
         cc=Gnew.faces.neighbors(cf,:);
         cc=cc(:);
         cc=unique(cc);
         if(any(cc==0))
            stopp=true;
         else
            [z,k]=min(Gnew.cells.z(cc));           
            if(z>=Gnew.cells.z(cell))
               stopp=true;
            else
               cell_line=[cell_line,cc(k)];
               cell_lines{count}=[cell_lines{count},cc(k)];
               cell=cc(k);
            end
         end
      end
      count=count+1;
      trap=traps(cell);
      J=[J;trap];
      if(trap==0)
         tstop=true;
      end
   end
   leaf_lines{leaf_count}=cell_line;
   leaf_count=leaf_count+1;
end
%
I=I(J>0);
J=J(J>0);
trap_matrix=sparse(I,J,1,max(traps),max(traps));