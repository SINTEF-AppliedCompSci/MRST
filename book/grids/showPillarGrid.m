grdecl = simpleGrdecl([4, 2, 3], .12, 'flat', true);
[X,Y,Z] = buildCornerPtPillars(grdecl,'Scale',true);
[x,y,z] = buildCornerPtNodes(grdecl);

%% Plot pillars
plot3(X',Y',Z','k'); drawAxisCross(.15); axis tight;
set(gca,'zdir','reverse'), view(35,35), axis off, zoom(1.2);
set(gca,'dataaspectRatio',[1.8,1.8,1])
print -deps2 showCPgrid-pillars.eps;

%% Plot points on pillars, mark pillars with faults red
hold on; I=[3 8 13];
hpr = plot3(X(I,:)',Y(I,:)',Z(I,:)','r','LineWidth',2);
hpt = plot3(x(:),y(:),z(:),'o'); 
hold off;
print -depsc2 showCPgrid-pts.eps;

%% Create grid and plot two stacks of cells
G = processGRDECL(grdecl);
args = {'FaceColor'; 'r'; 'EdgeColor'; 'k'};
hcst = plotGrid(G,[1:8:24 7:8:24],'FaceAlpha', .1, args{:});
print -dpng showCPgrid-cells.png;

%% Plot cells and fault surface
delete([hpt; hpr; hcst]);
plotGrid(G,'FaceAlpha', .15, args{:});
plotFaces(G, G.faces.tag>0,'FaceColor','b','FaceAlpha',.4);
print -dpng showCPgrid-grid.png;
