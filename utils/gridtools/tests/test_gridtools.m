%% Test data
G = cartGrid([10, 5, 3]);
plotPoints = @(pts) plot3(pts(:, 1), pts(:, 2), pts(:, 3), '.');

%% Test helpers
c = [1,2, 9, 7, 13, 51, 14, 40, 93];

% Test cells -> nodes function
[nodes, nodemap] = gridCellNodes(G, c);
coords = G.nodes.coords(nodes, :);

% Test cells -> faces function
[faces, facemap] = gridCellFaces(G, c);

%% Plot results
% Test selection of correct subset. By plotting faces along with the
% corresponding cells in different colors with alpha, we can see that they
% should overlap and form orange cells.
figure(1)
clf;
plotPoints(coords)
plotGrid(G, 'facec', 'none', 'edgea', .1);
plotGrid(G, c, 'facea', .5 , 'facec', 'yellow');
plotFaces(G, faces, 'facec', 'red');

% Test indexing into subset
figure(2)
clf;
target = 5;

ni = nodemap(target):nodemap(target+1)-1;
fi = facemap(target):facemap(target+1)-1;

assert(all(faces(fi) == indirectionSub(target, facemap, faces)));


plotPoints(coords(ni,:))
plotGrid(G, c(target), 'facea', .5);
plotGrid(G, 'facec', 'none', 'edgea', .1);
plotFaces(G, faces(fi), 'facec', 'red');
%% Plot coordinates corresponding to a subset of faces
figure(3)
clf;
f = [1, 3, 8,9, 23, 103, 49, 193];
[fnodes, fnmap] = gridFaceNodes(G, f);
plotPoints(G.nodes.coords(fnodes, :))
plotFaces(G, f);
plotGrid(G, 'facec', 'none', 'edgea', .1);

%% Some basic consistency checking on cellNo
cellno = gridCellNo(G);

fp = G.cells.facePos;
tmp = [];
for i = 1:numel(fp)-1
    v = unique(cellno(fp(i):(fp(i+1)-1)));
    % Each half face should belong to the same cell
    assert(numel(v) == 1);
    % All cells should have only one set of half faces
    assert(isempty(intersect(v, tmp)));
    tmp = [tmp v]; %#ok
end

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
