function [G, c] = producerInjectorGrid()
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

n = 30;
xMax = 1;
yMax = xMax;
mlqtMaxLevel   = 4;
wellGridFactor = 2^(-mlqtMaxLevel);


bnd = [0,0;0,xMax; xMax, xMax; xMax, 0];
d = xMax/n;
l = sqrt(2*(sqrt(2*xMax^2)/4)^2);
l1 = @(x) l - x;
l2 = @(x) 2*l - x;
l3 = @(x) 3*l - x;
% f = @(x,p1,p2) (x-p1(1))*(p2(2)-p2(1))/(p2(1)-p1(1)) + p1(2);
% f1 = @(x) f(x, [.1*xMax, .5*xMax], [.2*xMax, .9*xMax]);
offSet = 0.3*d;
fl = {[offSet l1(offSet); l1(offSet) offSet], ...
      [offSet l2(offSet); l2(offSet) offSet], ...
      [xMax-offSet, l3(xMax-offSet); l3(xMax-offSet), xMax-offSet]};
wl = {0.1*xMax*[1,1], 0.9*xMax*[1,1]};

pdims = [xMax, yMax];
%% Create sites
% Create well sites
wellGridSize = d * wellGridFactor;
[wellPts, wGs,protPts,pGs] = createWellGridPoints(wl, wellGridSize);

% Create fault sites
faultGridSize = d;
F = createFaultGridPoints(fl, faultGridSize);

% Create reservoir grid
dx = pdims(1)/ceil(pdims(1)/d);
dy = pdims(2)/ceil(pdims(2)/d);
vx = 0:dx:pdims(1);
vy = 0:dy:pdims(2);

[X, Y] = meshgrid(vx, vy);

resPtsInit = [X(:), Y(:)];

% Refine reservoir grid
if ~isempty(wellPts)
    res = {};
    varArg = {'level', 1, 'maxLev', mlqtMaxLevel, 'distTol', -1};
    for i = 1:size(resPtsInit,1)
        res = [res; mlqt(resPtsInit(i,:), wellPts, [d,d], varArg{:})];
    end
    resPts = vertcat(res{:,1});
    %resGridSize = 0.5*[res{:,2}]';
else
    resPts = resPtsInit;
    % resGridSize = repmat(0.5*min(dx,dy),size(resPts,1),1);
end

% Remove Conflic Points
resPts = removeConflictPoints2(resPts, wellPts,  wGs);
resPts = removeConflictPoints2(resPts, protPts,  pGs);
resPts = removeConflictPoints2(resPts, F.f.pts, F.f.Gs);
resPts = removeConflictPoints2(resPts, F.c.CC, F.c.R);

% Create Grid
Pts = [F.f.pts; wellPts;protPts; resPts];
G = clippedPebi2D(Pts, bnd);

% Tag grid
G.cells.tag = false(G.cells.num,1);
G.cells.tag(size(F.f.pts,1)+1:size(F.f.pts,1)+size(wellPts,1)) = true;

%%
G = computeGeometry(G);
c1 = G.cells.centroids(:,2) < l1(G.cells.centroids(:,1));
c2 = G.cells.centroids(:,2) < l2(G.cells.centroids(:,1)) & ~c1;
c3 = G.cells.centroids(:,2) < l3(G.cells.centroids(:,1)) & ~c1 & ~c2;
c4 = ~c1 & ~c2 & ~c3;
c = {c1, c2, c3, c4};

G.nodes.coords = 1000*G.nodes.coords;

%% Plot
% plotGrid(G)
% plotGrid(G,G.cells.tag,'facecolor','r')
%% save
% save producerInjectorGrid.mat

end
