function [z_spill_loc,g_trap,trap_level,z_spill_level,z_spill_loc_level,Gtop]=findTrappingStructure_dfs(Gt)
%addpath('/home/hnil/MATLAB/matlab_bgl')
addpath('/home/oddan/sim/matlab_bgl')
internal=all(Gt.faces.neighbors>0,2);
z_diff = Gt.cells.z(Gt.faces.neighbors(internal,2))-Gt.cells.z(Gt.faces.neighbors(internal,1));
top=accumarray([Gt.faces.neighbors(internal,2);Gt.faces.neighbors(internal,1)],double([z_diff;-z_diff]<0),[Gt.cells.num,1],@prod,0);
%A=sparse(double(Gt.faces.neighbors(internal,2)),double(Gt.faces.neighbors(internal,1)),z_diff>0,Gt.cells.num,Gt.cells.num);
%A=A+sparse(double(Gt.faces.neighbors(internal,1)),double(Gt.faces.neighbors(internal,2)),z_diff<0,Gt.cells.num,Gt.cells.num);
%A=speye(size(A))+A;
%clf,Gt.griddim=3;plotGrid(Gt,'FaceColor','none'),plotGrid(Gt,find(top),'FaceColor','r'),view(3)

top_cells=find(top);

Gtop=Gt;
 
 %mycells=sparse(top_cells,1:numel(top_cells),1,Gtop.cells.num,numel(top_cells));
 %{
 mycells=zeros(Gtop.cells.num,numel(top_cells));
 for i=1:numel(top_cells)
    mycells(top_cells(i),i)=1;
 end
 %}
%for kk=1:4
some_trapped=true;
trap_level={};
kk=0;
while some_trapped
   kk=kk+1;
   z_diff = Gtop.cells.z(Gtop.faces.neighbors(internal,2))-Gtop.cells.z(Gtop.faces.neighbors(internal,1));
   A=sparse(double(Gtop.faces.neighbors(internal,2)),double(Gtop.faces.neighbors(internal,1)),z_diff>=0,Gtop.cells.num,Gtop.cells.num);
   A=A+sparse(double(Gtop.faces.neighbors(internal,1)),double(Gtop.faces.neighbors(internal,2)),z_diff<=0,Gtop.cells.num,Gtop.cells.num);
   z_spill_level{kk}=nan(numel(top_cells),1);
   numcells_level{kk}=nan(numel(top_cells),1);
   trap_ind=[];
   
   for i=1:numel(top_cells)
       %disp('Test');    
      if(false) 
        seed=int32(top_cells(i));
        ct=int32(zeros(Gtop.cells.num,1));
        mc=int32(2);
        cc=dfs_jrn(A,seed-1,ct,mc);
        %tic;cc=dfs_jrn(A,seed-1,ct,mc);toc;
        m=find(cc>0);
        %clear cc;
      else
        cc=dfs(A',top_cells(i));
        cc=cc+1;cc(cc>0)=2;
        seed=top_cells(i);       
        m=find(cc>0);
         
         %clear cc
      end
      %{
      A_tmp=A';
      %tic;cc_tmp=bfs(A_tmp,top_cells(i));toc;
      %tic;cc_tmp=bfs_mex(A,top_cells(i),0);toc;
      tic;cc_tmp=dfs_mex(A,top_cells(i),0,0);toc;
      cc_tmp=cc_tmp+1;cc_tmp(cc_tmp>0)=2;      
      assert(all(cc==cc_tmp));
      %}
      Gtop.griddim=2;
      g=extractSubgrid(Gtop,m);
      %[g, gc, gf, gn] = extractSubgrid(Gtop, cells)
      
      k=any(g.faces.neighbors==0,2);
      %gb_cells=G.faces.neigbours(gf(k),
      cell=sum(g.faces.neighbors(k,:),2);
      [z_spill,j]=min(Gtop.cells.z(m(cell)));
      myloc_cc=(cc>0 & Gtop.cells.z<z_spill);
      if(any(myloc_cc))
         myloc_cc=(cc>0 & Gtop.cells.z<=z_spill);
         %myloc_cc(~internal)=0;
      end
      trap_ind=[trap_ind;find(myloc_cc),repmat(i,sum(myloc_cc),1)];
      tmp_mod= (cc>0 & Gtop.cells.z<=z_spill);
      Gtop.cells.z(find(tmp_mod))=z_spill;
      assert(sum(tmp_mod)>0);
      z_spill_level{kk}(i)=z_spill;
      %figure(2)
      %clf,pcolor(reshape(double(tmp_mod>0),Gtop.cartDims))
      %pause
  end
  disp(['Trap_level ',num2str(kk)])
   if(numel(top_cells)>0)
    trap_level{kk}=sparse(trap_ind(:,1),trap_ind(:,2),1,Gtop.cells.num,numel(top_cells));
    some_trapped=sum(trap_level{kk}(:))>0;
    %pcolor(reshape(full(double(sum(trap_level{kk},2)>0)),Gtop.cartDims));shading flat;colorbar
   else
     trap_level{kk}=sparse(Gtop.cells.num,numel(top_cells));
     some_trapped=false;
   end
end
z_spill_loc=zeros(Gtop.cells.num,1);%assuming zeros is above the top surface
for kk=1:numel(trap_level)
   z_spill_loc_level{kk}=z_spill_loc;%zeros(Gtop.cells.num,1);
   for i=1:numel(top_cells)
      ind=find(trap_level{kk}(:,i)>0);
      if(~isempty(ind))
         %z_spill_loc_level{kk}(ind)=max(z_spill_loc_level{kk}(ind),z_spill
         %_level{kk}(i));
         z_spill_loc_level{kk}(ind)=max(z_spill_loc(ind),z_spill_level{kk}(i));
         z_spill_loc(ind)=max(z_spill_loc(ind),z_spill_level{kk}(i));
      end
   end         
end
   
g_trap=trap_level{1};
for kk=2:numel(trap_level)
   %errind=sum(trap_level{kk}(trap_level{kk-1}>0) ~= 0,1);
   for i=1:size(trap_level{kk},2)
      ok=trap_level{kk}(trap_level{kk-1}(:,i)>0,i)==1;
      if(~all(ok))
         assert(all(trap_level{kk}(:,i)==0));
      end
   end
   %assert(all(trap_level{kk}(g_trap_level>0) ~= 0))
   %g_trap(g_trap==0 & trap_level{kk})=kk;
   g_trap(g_trap==0 &trap_level{kk})=kk;
end
%figure(33),clf,plotCellData(Gtop,z_spill_loc,find(z_spill_loc>0),'EdgeColor','none');plotGrid(Gtop,'FaceColor','none','EdgeAlpha',0.1)

%%
return
end
