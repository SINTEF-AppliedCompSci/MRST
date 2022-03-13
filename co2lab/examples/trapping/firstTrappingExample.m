%% Trapping example
% In this example, we demonstrate the basic routines for computing traps,
% accumulation areas, spill paths, etc. To this end, we will use a simple
% shoe-box with a dip and a perturbed top surface.

mrstModule add co2lab coarsegrid
%% Make grid and rock data
% We consider a sandbox of dimensions 10 km x 5 km x 50 m, that lies at a
% depth of 1000 m and has an inclination in the x-direction. For the top
% surface, we add a smooth sin/cos perturbation that will create domes. The
% porosity is set uniformly to 0.25
[Lx,Ly,H] = deal(10000, 5000, 50);
G  = cartGrid([100 100 1],[Lx Ly H]);
x  = G.nodes.coords(1:G.nodes.num/2,1)/Lx; 
y  = G.nodes.coords(1:G.nodes.num/2,2)/Ly;
z  = G.nodes.coords(1:G.nodes.num/2,3)/H;
zt = z + x - 0.2*sin(5*pi*x).*sin(5*pi*y.^1.5) - 0.075*sin(1.25*pi*y) + 0.15*sin(x+y);
zb = 1 + x;
G.nodes.coords(:,3) = [zt; zb]*H+1000;
G = computeGeometry(G);

%% Construct top surface grid
% We create a hybrid grid that represents the top surface. The grid is a 2D
% grid defined as a surface in 3D and contains a mapping from the 2D cells
% on the surfrace to the column of volumetric cells that lie below in the
% 3D grid. For visualization purposes, we create an extra grid that lies
% above the true top surface.
Gt = topSurfaceGrid(G);

figure;
Gt_zshifted = Gt; 
Gt_zshifted.nodes.z = Gt_zshifted.nodes.z - 15;
plot_opts = {'edgeColor', 'k', 'edgeAlpha', 0.1};
plotGrid(G, plot_opts{:});
plotCellData(Gt_zshifted, Gt_zshifted.cells.z, plot_opts{:});
view(30,25); axis tight

%% Geometric analysis of caprock (spill-point analysis)
% Compute traps and spill paths connecting them. Here, we use the
% cell-based algorithm. The cells that belong to the identified traps are
% colored white in the plot
res = trapAnalysis(Gt, true);
num_traps = max(res.traps);
plotGrid(Gt_zshifted, res.traps>0, 'FaceColor','white', 'EdgeAlpha',.1);


%% Show connection between traps and spill paths
% We make a 2D plot of the top surface in which traps are colored red,
% cells that lie along the connecting spill paths are colored green, and
% the remaining cells are colored blue. In addition, we display the number
% associated with each trap slightly above its local maximum.
clf
fpos = get(gcf,'Position');
set(gcf,'Position',[300 400 800 420],'PaperPositionMode','auto');
trap_field = zeros(size(res.traps));
trap_field(res.traps>0) = 2;
for r = [res.cell_lines{:}]'
    for c = 1:numel(r)
        trap_field(r{c}) = 1;
    end
end

subplot(2,3,[1 4]);
plotCellData(Gt, trap_field, 'EdgeColor','none');
view(90,90); axis tight off

for i=1:num_traps
   ind = res.top(i);
   text(Gt.cells.centroids(ind,1),Gt.cells.centroids(ind,2), res.trap_z(i)-50,...
      num2str(res.traps(ind)), 'Color',.99*[1 1 1], ...
      'FontSize', 14, 'HorizontalAlignment','center');
end
colormap(jet);

%% Compute the total trapping capacity
% We compute the pore volume inside each of the traps, and report in
% descending order. In addition, we relate the trap volume to the total
% volume of the sandbox.
rock2D.poro = 0.25 * ones(G.cells.num, 1);
trap_volumes = volumesOfTraps(Gt, res, 1:num_traps, 'poro', rock2D.poro);

total_trapping_capacity = sum(trap_volumes);
pv = sum(poreVolume(Gt,rock2D));
fprintf('Total trapping capacity is: %6.3e\n', total_trapping_capacity);
fprintf('This amounts to %.2f %% of a total pore volume of %6.2e\n\n', ...
   total_trapping_capacity/pv*100, pv);
[sorted_vols, sorted_ix] = sort(trap_volumes, 'descend');

fprintf('trap ix   | trap vol(m3)  | cells in trap\n');
fprintf('----------+---------------+--------------\n');
tcells = zeros(num_traps,1);
for i = 1:num_traps
   tcells(i) = sum(res.traps == sorted_ix(i));
   fprintf('%7d   |  %10.3e   | %5d\n', ...
       sorted_ix(i), sorted_vols(i), tcells(i));
end
fprintf(['\nTogether, the five largest traps cover %6.2e m3, which represents %3.1f%% of' ...
   '\nthe total trapping capacity of this grid.\n'], ...
   sum(sorted_vols(1:5)), sum(sorted_vols(1:5)) / total_trapping_capacity * 100)

subplot(2,3,[2 3]);
bar(1:num_traps,sorted_vols);
set(gca,'XTickLabel',sorted_ix);
title('Trap volume (m^3)');

subplot(2,3,[5 6]);
bar(1:num_traps,tcells);
set(gca,'XTickLabel',sorted_ix);
title('Number of cells in trap');

%% Plot leaf nodes
% Next, we will look at leaf nodes in our trapping system, i.e., traps that
% only have outgoing 'rivers'. The leaf traps are shown in white color. The
% drainage areas and the resulting migration path are given a unique  color
% so that one can trace the migration of CO2 upward until it either  leaves
% the formation or reaches a higher trap that accumulates CO2 from multiple
% leaf traps. All traps that are not leaf traps are shown in dark gray.
% To do this, we will utilize two functions that are called as part of
% trapAnalysis, but whose retun argument are not fully exposed by the
% higher-level interface in trapAnalysis.
moduleCheck('matlab_bgl');
ts         = findTrappingStructure(Gt);
trap_con   = findTrapConnections(Gt, ts.z_spill_loc);
traps      = trap_con.traps;
leaf_lines = trap_con.leaf_lines;
leaf_traps = trap_con.leaf_traps;

% Plot traps and grid
clf, set(gcf,'Position',fpos);
plotGrid(ts.Gtop, ts.z_spill_loc>0, 'EdgeColor', 'none', 'FaceColor',[.25 .25 .25])
plotGrid(ts.Gtop, 'FaceColor', 'none', 'EdgeAlpha', .1)

% Find leaves and the corresponding traps
[ll, tc, dc] = deal([]);
C    = maxTPFAGravityMatrix(ts.Gtop);
nll  = numel(leaf_lines);
p    = randperm(nll);
lvol = zeros(nll,1);
tic
for k=nll:-1:1
   % extract the leaf line out of the trap
   ll  = [ll; leaf_lines{k}', p(k)*ones(length(leaf_lines{k}),1)]; %#ok
   
   % find all cells in the given trap
   t_cells = find(traps ==leaf_traps{k});
   tc = [tc; t_cells]; %#ok

   % find all cells in the drainage region of the leaf
   cc = find( dfs(C', double(t_cells(1))) >0);
   dc  = [dc; cc, p(k)*ones(length(cc),1)]; %#ok
   
   % calculate the leaf volume
   cc = find( dfs(C, double(t_cells(1)))>0 );
   lvol(k) = sum((ts.Gtop.cells.z(cc) - Gt.cells.z(cc)).*Gt.cells.volumes(cc));
   disp(['Leaf volume : ', num2str(lvol(k))]);
end

% Plot leaves and the deepest traps
plotCellData(ts.Gtop, [ll(:,2); dc(:,2)], [ll(:,1); dc(:,1)], 'EdgeColor', 'none');
plotGrid(ts.Gtop, tc, 'EdgeColor', 'none', 'FaceColor', [.9 .9 .9]);
set(gca,'Projection','perspective');
view(120,25); axis tight


%% Plot accumulation areas using interactive viewer
% Finally, we launch the interactive viewer. Here, each trap and the
% surrounding accumulation areas are given a unique color. Injection
% scenarios for a single well can be investigated by clicking on the
% surface (inside a trap or an accumulation region). This will compute the
% resulting spill-point path and the potential trapping volumes
% inside traps that are visited before the spill path exits the model.
h = interactiveTrapping(Gt, 'method', 'cell', 'light', true, ...
   'spillregions', true, 'colorpath', false, 'injpt', 1898);
view(80,20);
set(gca,'DataAspect', [2 1.5 .02]);

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
