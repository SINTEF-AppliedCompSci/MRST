mrstModule add upr dg vem vemmech matlab_bgl agglom coarsegrid incomp

%% Common params

% Figures 
pos = [-1000, 0, 500, 500];
fontSize = 12;
pth = fullfile(mrstPath('dg'), 'examples', 'rsc-2019', 'cubature', 'fig');
if 1
    savepng = @(name) print(fullfile(pth, name), '-dpng', '-r300');
    saveeps = @(name) print(fullfile(pth, name), '-depsc');
else
    savepng = @(name) [];
    saveeps = @(name) [];
end

cNo = 7;
gray = [1,1,1]*0.8;
cw = @(w,wMax,wMin) 20/(wMax - wMin)*(w - wMin) + 15;
clr = lines(2);

%% make PEBI grid

n  = 10;
l = 1000;
G = pebiGrid(l/n, [l,l]);

t = pi/12;
R = [cos(t), -sin(t); sin(t), cos(t)];
G.nodes.coords = G.nodes.coords*R';

close all
plotGrid(G)
axis equal tight
G = computeGeometry(G);
G = computeVEMGeometry(G);
G = computeCellDimensions2(G);
[G.cells.equal, G.faces.equal] = deal(false);

%% Make coarse grid

rock = makeRock(G, 1,1);
T = computeTrans(G, rock);
p = partitionMETIS(G, T, 15);
GC = generateCoarseGrid(G, p);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC, 'useFullBB', true);
[GC.cells.equal, GC.faces.equal] = deal(false);

%%

internalConn = ~any(GC.faces.neighbors == 0,2);
cubMom  = MomentFitting2DCubature(GC, 3, internalConn);
cubFull = MomentFitting2DCubature(GC, 3, internalConn, 'reduce', false);
cubCoarse = CoarseGrid2DCubature(GC, 3, internalConn);
cubTri = TriangleCubature(G, 3, []);

%%

[~, x, w1] = cubCoarse.getCubature(cNo, 'volume');
[~, x, w2] = cubMom.getCubature(cNo, 'volume');
wMax = max([w1; w2]);
wMin = min([w1; w2]);
cw = @(w) 50/(wMax - wMin)*(w - wMin) + 15;

%%

close all
figure('Position', pos);
hold on
plotGrid(GC, cNo, 'facec', gray);
c = GC.partition == cNo;
plotGrid(G, c, 'linestyle', '--', 'facec', 'none');
box on;
axis equal tight
ax = gca;
[ax.XTick, ax.YTick] = deal([]);

hold on
faces = GC.cells.faces(GC.cells.facePos(cNo):GC.cells.facePos(cNo+1)-1);
for f = faces
    xf = GC.nodes.coords(GC.faces.nodes(GC.faces.nodePos(f):GC.faces.nodePos(f+1)-1),:);
    plot(xf(:,1), xf(:,2), '-', 'color', clr(1,:), 'linew', 2);
end
hold off
saveeps('coarse-details');

%%

close all
figure('Position', pos);
hold on
plotGrid(GC, cNo, 'facec', gray);
box on;
axis equal tight
ax = gca;
[ax.XTick, ax.YTick] = deal([]);
saveeps('coarse');

%%

figure('Position', pos);
hold on
plotGrid(GC, cNo, 'facec', gray);
c = find(GC.partition == cNo);

triNo = mcolon(cubTri.triangulation.triPos(c), cubTri.triangulation.triPos(c+1)-1);
for t = triNo
    ix = cubTri.triangulation.ConnectivityList(t,:);
    xt = cubTri.triangulation.Points(ix,:);
    xt = [xt; xt(1,:)];
    plot(xt(:,1), xt(:,2), 'k--');
end

% plotGrid(G, c, 'facec', 'none', 'linestyle', 'k--');
[~, x, w] = cubCoarse.getCubature(cNo, 'volume');
wMax = max(w);
wMin = min(w);
for pNo = 1:numel(w)
    plot(x(pNo,1), x(pNo,2), '.', 'markerSize', cw(w(pNo)), 'color', 'k');
end


box on;
axis equal tight
ax = gca;
[ax.XTick, ax.YTick] = deal([]);
saveeps('coarse-tri');

%%

% close all
figure('Position', pos);
plotGrid(GC, cNo, 'facec', gray);
hold on
[~, x, w] = cubFull.getCubature(cNo, 'volume');

wMax = max(w);
wMin = min(w);
nClr = 40;
for pNo = 1:numel(w)
%     cIx = floor((w(pNo) - wMin)./(wMax - wMin)*(nClr-1)) + 1;
    plot(x(pNo,1), x(pNo,2), '.', 'markerSize', 30, 'color', 'k');
end
axis equal tight
box on;
axis equal tight
ax = gca;
[ax.XTick, ax.YTick] = deal([]);
saveeps('coarse-init')

%%

% close all
figure('Position', pos);
plotGrid(GC, cNo, 'facec', gray);
hold on
[~, x, w] = cubMom.getCubature(cNo, 'volume');

clr = jet(nClr);
wMax = max(w);
wMin = min(w);
nClr = 40;
for pNo = 1:numel(w)
    cIx = floor((w(pNo) - wMin)./(wMax - wMin)*(nClr-1)) + 1;
    plot(x(pNo,1), x(pNo,2), '.', 'markerSize', cw(w(pNo)), 'color', 'k');
end

axis equal tight
axis equal tight
box on;
axis equal tight
ax = gca;
[ax.XTick, ax.YTick] = deal([]);
saveeps('coarse-reduced')


%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
