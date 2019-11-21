clear
close all
%%
G = cartGrid([25, 25], [200, 200]);
G = computeGeometry(G);

ns = 10;
th = linspace(0.18, 0.32, ns+1)' * pi;
R = 150;
pW = [R*cos(th), R*sin(th)];

pbdy = [50,150; 60,80;120,36;160,60;150,110;100,170];

cCtro = G.cells.centroids;
in = inpolygon(cCtro(:,1), cCtro(:,2), pbdy(:,1), pbdy(:,2));
CI = find( in );
CO = find(~in );

[bn, bc] = demo_getBdyNodesCells(G, CI);

figure, hold on, axis equal tight off
plotGrid(G, CO, 'facecolor', 'none')
plotGrid(G, CI, 'facecolor', 'y')
demo_plotLine(pW,   'ko-', 'b', 4)
demo_plotPoly(pbdy, 'k^-', 'r', 5)
demo_plotPoly(G.nodes.coords(bn,:), 'ks-', 'g', 4)
legend('G - CO', 'G - CI', 'Well path' ,'Specified VOI boundary','Clipped VOI boundary')
%%
ly = 15;
ny = 8;
na = 6;
pWR0 = arrayfun(@(ii)pointsSingleWellNode(pW, ly, ny, na, ii), (1:ns+1)');
pWR  = [vertcat(pWR0.cart); vertcat(pWR0.rad)];
assert(all(inpolygon(pWR(:,1), pWR(:,2), pbdy(:,1), pbdy(:,2))),...
    ['Points outside the boundary were detected, please reduce ',...
    'the size of Cartesian region']);
[tWR, ~, bnWR] = getConnListAndBdyNodeWR2D(pWR0, ny, na);

figure, hold on, axis equal tight off
gWR = tessellationGrid(pWR, tWR);
plotGrid(gWR, 'facecolor', 'none')
plotGrid(G, CO, 'facecolor', 'none')
demo_plotPoly(G.nodes.coords(bn, :), 'ks-', 'g', 4)
demo_plotPoly(pWR(bnWR, :), 'ko-', 'y', 3)
legend('G', 'gWR (WR grid)', 'VOI boundary' ,'WR boundary')

