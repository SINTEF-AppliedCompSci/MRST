%% Inspect the Hugin West Formation
% This formation has the largest difference in the trap volumes. We extract
% a subset in the south of the model and compare the traps computed by the
% node-based and the cell-based algorithms
mrstModule add libgeometry deckformat coarsegrid
grdecl = getAtlasGrid('Huginfmwest');
try
   G = processgrid(grdecl{1}); clear grdecl;
   G = mcomputeGeometry(G(1));
catch
   G = processGRDECL(grdecl{1}); clear grdecl;
   G = computeGeometry(G(1));
end


%% Show the whole formation with trapping structure
Gt = topSurfaceGrid(G);
interactiveTrapping(Gt, 'method', 'node', 'light', true, 'spillregions', true);
set(gca,'Position',[.075 .075 .85 .85])
view(-70,55)
axis on
set(gca,'XMinorGrid', 'on', 'YMinorGrid', 'on', 'ZMinorGrid', 'on', ...
   'YTickLabel',[],'XTickLabel',[]);

%% Extract a subset
cDims = G.cartDims;
cntrd = G.cells.centroids;
G     = extractSubgrid(G, find((cntrd(:,1) > 4.35e5) & (cntrd(:,2)<6.45e6)));
G.cartDims = cDims;
Gt = topSurfaceGrid(G);
clear cntrd cDims;

% Initialize plotting
clf;
fcol   = get(gcf,'Color');
map    = (1*hsv(6) + 1.5*repmat(fcol, 6, 1))./2.5;
map    = map([1 1:end],:); map(1,:) = fcol;
p      = get(gcf,'Position'); set(gcf,'Position', [p(1:2) 840 420]);

%% Node-based traps
% Find traps using the node-based algorithms. To get a consistent ordering
% between the two algorithms, we use a perturbation vector to renumber the
% traps, from back to front in the plot.
ta  = trapAnalysis(Gt, false);
pn  = [0 5 2 3 6 1 4]';
[~, pn_inv] = sort(pn);
val = ta.trap_regions;

% Plot accumulation regions
subplot(1,2,1); cla
h   = plotCellData(Gt, ones(Gt.cells.num,1),'EdgeColor', 'none');
%set(h, 'FaceVertexCData', map(pn(val+1)+1,:));
set(h, 'FaceVertexCData', map(pn_inv(val+1),:));

% Plot traps
%plotCellData(Gt, pn(ta.traps+1), (ta.traps ~= 0), 'EdgeColor', 'k')
plotCellData(Gt, pn_inv(ta.traps+1)-1, (ta.traps ~= 0), 'EdgeColor', 'k')

% Fix axis, set light and colorbar
view(-20,40); axis tight
light('Position',[-1 -1 -1],'Style','infinite');lighting phong
colormap(hsv(6));
colorbar('horiz'); caxis([.5 6.5]);
set(gca,'XTickLabel',[],'YTickLabel',[]);

% Make table with the number of traps and corresponding volumes as computed
% by the node-based algorithm
tcn = accumarray(ta.traps+1,1);
tvn = volumesOfTraps(Gt,ta,1:6);
fprintf('\n\nNode-based method:\n');
fprintf('Trap     ');  fprintf('& %7d ', 1:6); fprintf('\\\\\\hline\n');
fprintf('Cells    '); fprintf('& %7d ', tcn(pn(2:end)+1)); fprintf('\\\\\n');
fprintf('Volume   '); fprintf('& %4.1e ',tvn(pn(2:end))); fprintf('\\\\\\hline\n\n');

%% Cell-based traps
% Find traps using the node-based algorithms. To get a consistent ordering
% between the two algorithms, we use a perturbation vector to renumber the
% traps, from back to front in the plot.
ta  = trapAnalysis(Gt, true);
%pc  = [0 1 4 2 5 3 6]';
pc  = [0 1 3 5 2 4 6]';
[~, pc_inv] = sort(pc);
val = ta.trap_regions;

% Plot accumulation regions
subplot(1,2,2); cla
h   = plotCellData(Gt, ones(Gt.cells.num,1), 'EdgeColor', 'none');
%set(h, 'FaceVertexCData', map(pc(val+1)+1,:));
set(h, 'FaceVertexCData', map(pc_inv(val+1),:));

% Plot traps
%plotCellData(Gt, pc(ta.traps+1), (ta.traps ~= 0), 'EdgeColor', 'k')
plotCellData(Gt, pc_inv(ta.traps+1)-1, (ta.traps ~= 0), 'EdgeColor', 'k')

% Fix axis, set light and colorbar
view(-20,40); axis tight
light('Position',[-1 -1 -1],'Style','infinite');lighting phong
colorbar('horiz'); caxis([.5 6.5]);
set(gca,'XTickLabel',[],'YTickLabel',[]);

% Make table with the number of traps and corresponding volumes as computed
% by the cell-based algorithm
tcc = accumarray(ta.traps+1,1);
tvc = volumesOfTraps(Gt,ta,1:6);
fprintf('\n\nCell-based method:\n');
fprintf('Trap     '); fprintf('& %7d ', 1:6); fprintf('\\\\\\hline\n');
fprintf('Cells    '); fprintf('& %7d ', tcc(pc(2:end)+1)); fprintf('\\\\\n');
fprintf('Volume   '); fprintf('& %4.1e ',tvc(pc(2:end))); fprintf('\\\\\\hline\n\n');

%%
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
