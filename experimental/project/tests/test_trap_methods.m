coarsening=5
[utsira_grdecl,utsira_meta]=getAtlasGrid('Utsirafm','coarsening',coarsening)
utsira_grdecl=utsira_grdecl{1};
G3D_utsira=processGRDECL(utsira_grdecl);
p=0;num_cells=0;
for i=1:numel(G3D_utsira)
    if(G3D_utsira(i).cells.num>num_cells)
        num_cells=G3D_utsira(i).cells.num;
        p=i;
    end
    
end
g_utsira=topSurfaceGrid(G3D_utsira(p));
g_utsira.cells.cellNodes=getCellNodes(g_utsira);
profile off;
profile on;
use_multipoint=false;
trap_str=findTrappingStructureMaxUpwind(g_utsira,'use_multipoint',use_multipoint);
%trap_str=findTrappingStructure(g_utsira,'use_9p',true);
z_spill_loc1=trap_str.z_spill_loc;
trap_str=findTrappingStructure(g_utsira,'use_multipoint',use_multipoint);
z_spill_loc=trap_str.z_spill_loc;
use_cell_based=false
trap_str_new = trapAnalysis(g_utsira, use_cell_based);
tmp=[0;trap_str_new.trap_z];
z_spill_loc_new=tmp(trap_str_new.traps+1);
profile off;
profile viewer
[sommet, region] = spill_field_cells(g_utsira);
Gtop=trap_str.Gtop;
%%
g_utsira.cells.cellNodes=getCellNodes(g_utsira);
figure(1),clf
subplot(1,3,1)
plotGrid(g_utsira,z_spill_loc1>0,'EdgeColor','b')
plotGrid(g_utsira,'EdgeColor','k','FaceColor','none')
subplot(1,3,2)
plotGrid(g_utsira,z_spill_loc>0,'EdgeColor','r')
plotGrid(g_utsira,'EdgeColor','k','FaceColor','none')
subplot(1,3,3)
plotGrid(g_utsira,z_spill_loc_new>0,'EdgeColor','r')
plotGrid(g_utsira,'EdgeColor','k','FaceColor','none')
%%
figure(2),clf
val=double(z_spill_loc>0)+double(z_spill_loc1>0)*2+double(z_spill_loc_new>0)*2^2;
%{
plotGrid(g_utsira,z_spill_loc>0,'EdgeColor','none')
plotGrid(g_utsira,z_spill_loc1>0,'EdgeColor','none','FaceColor','b')
plotGrid(g_utsira,(z_spill_loc1>0)~=(z_spill_loc>0),'EdgeColor','none','FaceColor','r')
plotGrid(g_utsira,'FaceColor','none','edgea',0.1)
%}
plotCellData(g_utsira,val,find(val>0),'edgea',0.1)
plotGrid(g_utsira,'FaceColor','none','edgea',0.1)
colorbar
axis tight;axis off;view(2);%view([1 -4 6])
drawnow;
