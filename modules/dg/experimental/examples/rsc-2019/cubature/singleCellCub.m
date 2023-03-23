mrstModule add upr dg vem vemmech

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

cNo = 100;
gray = [1,1,1]*0.8;
cw = @(w,wMax,wMin) 20/(wMax - wMin)*(w - wMin) + 15;

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

%%

internalConn = ~any(G.faces.neighbors == 0,2);
cubMom  = MomentFitting2DCubature(G, 3, internalConn);
cubFull = MomentFitting2DCubature(G, 3, internalConn, 'reduce', false);
cubTri  = TriangleCubature(G, 3, internalConn);

%%

[~, x, w1] = cubTri.getCubature(cNo, 'volume');
[~, x, w2] = cubMom.getCubature(cNo, 'volume');
wMax = max([w1; w2]);
wMin = min([w1; w2]);
cw = @(w) 50/(wMax - wMin)*(w - wMin) + 15;

%%

close all
figure('Position', pos);
hold on
plotGrid(G, cNo, 'facec', gray);
box on;
axis equal tight
ax = gca;
[ax.XTick, ax.YTick] = deal([]);
saveeps('pebi');

%%

figure('Position', pos);
hold on
plotGrid(G, cNo, 'facec', gray);
[~, x, w] = cubTri.getCubature(cNo, 'volume');
wMax = max(w);
wMin = min(w);
for pNo = 1:numel(w)
    plot(x(pNo,1), x(pNo,2), '.', 'markerSize', cw(w(pNo)), 'color', 'k');
end

f = G.cells.faces(G.cells.facePos(cNo):G.cells.facePos(cNo+1)-1);
n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
n = reshape(n,2,[])';
sgn = 1 - 2*(G.faces.neighbors(f,1) ~= cNo);
n(sgn<0,:) = n(sgn<0,[2,1]);
n = n(:,1);
xn = G.nodes.coords(n,:);
xc = G.cells.centroids(cNo,:);
for nNo = 1:numel(n)
    xx = [xn(nNo,:); xc];
    plot(xx(:,1), xx(:,2), '--k');
end

box on;
axis equal tight
ax = gca;
[ax.XTick, ax.YTick] = deal([]);
saveeps('pebi-tri');

%%

% close all
figure('Position', pos);
plotGrid(G, cNo, 'facec', gray);
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
saveeps('pebi-init')

%%

% close all
figure('Position', pos);
plotGrid(G, cNo, 'facec', gray);
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
saveeps('pebi-reduced')

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
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
