
% to do trap analysis we need  the following modules
%mrstModule add internal/vertical-equil
%mrstModule add gridprocessing
%mrstModule add mex/libgeometry
%mrstModule add mex
%mrstModule add /home/hnil/heim/SVN/simmatlab/branches/halvor/mrst_extra/gridprocessing
%addpath('/home/hnil/heim/SVN/simmatlab/branches/halvor/mrst_extra/gridprocessing')
mycase='igems_3d'
mycase='simple'
% mycase='igems_surface'
switch mycase
  case 'sara_3d'
    %%
    if(false)
        a=load('/home/hnil/heim/SVN/simmatlab/branches/mrst-reorg/modules/internal/vertical-equil/examples/igems/data/sara/g_top_FMM.mat');
        Gt=a.g_top_FMM;
    else
        a=load('/home/hnil/heim/SVN/simmatlab/branches/mrst-reorg/modules/internal/vertical-equil/examples/igems/data/sara/g_top_OSS.mat');
        Gt=a.g_top_OSS;          
    end
    Gt.griddim=2;
    G.cells.cellNodes=getCellNodes(Gt);
    clear a;
  case 'igems_3d'
    mygrid='FMM'
    %mymod='z1_coarse44_uniform'
    mymod='z1_uniform'
    [Gt,rock2D,SVE,rock,invert_axis,G]=readIGEMSgrid(mygrid,mymod);
    
  case 'igems_surface'
    % we can load surfaces from the igems data set the possible names are
    % the directories  in data/surfaces 
    name='flatNP2';
    name='OSSNP2';
    % there is 100 different surfaces we choose number
    i=2;
    % For testing it is usefull to use a coarse variant
    coarse=[1,1];
    %[Gt,rock2D,SVE,rock,G]=readIGEMSIRAP(name,i,coarse)
    [Gt,rock2D,rock,G]=readIGEMSIRAP(name,i,coarse)
  case 'simple'
    % simple test grid for finding traps
    L=1000;L_p=L/5;H=10
    G=cartGrid([101 101 1],[L L H]);
    G.nodes.coords(:,3)=100+G.nodes.coords(:,3)+G.nodes.coords(:,1)*0.01...
        -2*sin(pi*G.nodes.coords(:,1)/L_p).*sin(pi*G.nodes.coords(:,2)/L_p);
    G=computeGeometry(G);
    Gt=topSurfaceGrid(G);
  otherwise
    error()
end
%% find traps
tic;[z_spill_loc,g_trap,trap_level,z_spill_level,z_spill_loc_level,Gnew]=findTrappingStructure_dfs(Gt);toc;
disp('Finnished finding traplevels')
% find the connection between traps
[cell_lines,traps]=findCellLines(Gnew,z_spill_loc);
disp('End find rivers')
%% plot traps

if (false)
    % plot traps and connection on  the flat surface
    figure(1)
    clf,a=plotCellData(Gnew,traps,traps>0)
    plotGrid(Gnew,'FaceColor','none')
    plotGrid(Gnew,[cell_lines{:}],'FaceColor','k');
    set(gca,'XTick',[],'YTick',[]);
    axis equal tight
    axis off
else
    % plot traps and connection on  the flat surface
    figure(11)
    %clf,a=plotCellData(Gt,traps,traps>0)
    clf,a=plotCellData(Gt,z_spill_loc,z_spill_loc>0)
    plotGrid(Gt,'FaceColor','none','EdgeColor','none')
    plotGrid(Gt,[cell_lines{:}],'FaceColor','k');
    axis equal tight
    axis off
    view(3)
end
%%
%plotGrid(Gnew,'FaceColor','none')
%plotGrid(Gnew,[cell_lines{:}],'FaceColor','k');
set(gca,'XTick',[],'YTick',[]);
axis equal tight
axis off

%%
% plot cell lines
figure(2)
clf,plotCellData(Gt,Gt.cells.z,'EdgeAlpha',0.1,'EdgeColor','k')
view(3)
axis off
%%
set(gca,'XTick',[],'YTick',[],'ZTick',[])
return
figure(1),clf
pcolor(reshape(z_spill_loc,Gt.cartDims));shading interp
rivers=zeros(Gt.cells.num,1);
rivers([cell_lines{:}])=1;
% add lakes
rivers(z_spill_loc>0)=1;
figure(2),clf
pcolor(reshape(rivers,Gt.cartDims));shading interp

%hold on

