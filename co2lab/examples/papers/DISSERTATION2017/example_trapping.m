%% Loading the Utsira grid 
coarsening_level = 1; % grid downsampling factor (1 = no downsampling)
Gt = getFormationTopGrid('Utsirafm', coarsening_level);

% We check the number of cells in the top surface grid.
Gt.cells.num

%% Visualize grid

% Zenith view
plotCellData(Gt, Gt.cells.z, 'edgecolor', 'none')
view(90, 90); axis tight;
c = colorbar;  c.Label.String = 'Depth (m)';
set(c, 'fontsize', 20); set(gca, 'fontsize', 20);
set(gcf, 'position', [0 0 1400 600]); camlight headlight; material dull

% Inclined view
view(80, 50);

%% Compute trapping structure
% The second argument indicates node-based (as opposed to cell-based)
% interpretation of geometry
ta = trapAnalysis(Gt, false)

% show number of traps
max(ta.traps)

%% Visualize the traps

% Zenith view
clf;
plotCellData(Gt, ta.traps, 'edgecolor', 'none');
view(90, 90); axis tight; colormap lines;
set(gca, 'fontsize', 20);
set(gcf, 'position', [0 0 1400 600]); camlight headlight; material dull

% Inclined  view
view(80, 50); 

%% Visualize spill regions and rivers

% construct color field identifying traps
trapfield = ta.traps;
trapfield(trapfield==0) = NaN;

% construct color field identifying rivers
river_field = NaN(size(ta.traps));
for i = 1:numel(ta.cell_lines)
   rivers = ta.cell_lines{i};
   for r = rivers
      river_field(r{:}) = 1;
   end
end

% Top view

% plot traps in solid color
clf; plotCellData(Gt, trapfield, 'edgecolor', 'none'); 

% plot regions in semi-transparent color
plotCellData(Gt, ta.trap_regions, 'facealpha', 0.5, 'edgecolor', 'none');

% add rivers
plotCellData(Gt, river_field, 'edgealpha', 0);

view(90, 90); axis tight; colormap lines;
set(gca, 'fontsize', 20);
set(gcf, 'position', [0 0 1400 600]); camlight headlight; material dull

% Inclined view
view(80, 50);

%% Topographical map plot

h = figure;
mapPlot(h, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines);
set(gca, 'fontsize', 20);
set(gcf, 'position', [0 0 1400 600]); view(90, 90);

%% Trap capacity in volumetric and mass terms

% We add a zero to the vector of spill point depths.  This is necessary for
% indexing purposes, as explained in the next comment.
depths = [0; ta.trap_z];

% The line below gives the correct "trap height" value for cells within
% traps.  For the other cells, the value is incorrect, but will be fixed by
% the succeeding line.  Note that since we have added an entry to the front
% of the 'depths' vector, we increase indices by 1.  (Otherwise, the
% presence of zeros in 'ta.traps' would cause an indexing error, since Matlab
% arrays are indexed from 1, not 0).
h = min((depths(ta.traps+1) - Gt.cells.z), Gt.cells.H);
h(ta.traps==0) = 0; % cells outside any trap should have zero trap volume

% Bulk cell volumes are now simply obtained by multiplying by the respective
% areas.  For non-trap cells, the volumes will be zero, since 'h' is zero for
% these cells.
cell_tvols = h .* Gt.cells.volumes;

% We accumulate the trap volume of individual cells to obtain the total bulk
% trap volume for each structural tra
tvols = accumarray(ta.traps+1, cell_tvols);

% The first entry of 'tvols' contain the combined trap volume of all non-trap
% cells.  Obviously, this value should be zero, but let us check this to make
% sure.
assert(tvols(1)==0);

% We remove the first entry of 'tvols', since we are only interested in the
% volumes of the actual traps.
tvols = tvols(2:end);


porosity = 0.36; 
seafloor_temp = 7 + 273.15; % Seafloor temperature, in degrees Kelvin
seafloor_depth = 80; % in meters
temp_grad = 35.6; % temperature increase per kilometer depth (deg. Kelvin)
rho_brine = 1000; % brine density (kilogram per cubic meter)

% Ensure that gravity is not zero
gravity on;

% computing temperature field (at caprock level)
T = seafloor_temp + temp_grad .* (Gt.cells.z - seafloor_depth)/1000;

% computing hydrostatic pressure (at caprock level)
P = 1 * atm + rho_brine * norm(gravity) * Gt.cells.z;

% Making stripped-down fluid object that only contains CO₂ density function
fluid = addSampledFluidProperties(struct, 'G');

% Computing local CO₂ (at caprock level)
CO2_density = fluid.rhoG(P, T);

% Computing trap capacity in mass term for each cell
cell_tmass = cell_tvols .* porosity .* CO2_density;

% Accumulating cell values to get trap capacity in mass term for each trap
tmass = accumarray(ta.traps+1, cell_tmass);
assert(tmass(1)==0);

% We produce a bar plot where trap volumes are plotted in descending order.
% Here too, we notice the presence of a handful of large traps, and a long
% tail of vanishingly small traps.
clf;
subplot(1,2,1);
bar(sort(tvols, 'descend'))
xlabel('Traps (sorted by volume)');
ylabel('Cubic meters'); set(gca, 'fontsize', 20);
title('Capacity (volume)');
subplot(1,2,2);
bar(sort(tmass/1e9, 'descend'));
xlabel('Traps (sorted by mass capacity)');
ylabel('Megatonnes'); set(gca, 'fontsize', 20);
title('Capacity in mass terms');
set(gcf, 'position', [0 0 1300 440]);


% Plot traps according to capacity in mass terms (megatonnes)
figure
plotCellData(Gt, tmass(ta.traps+1) / 1e9, 'edgealpha', 0.1);
view(90, 90); axis tight; colormap cool; 
c = colorbar; c.Label.String = 'Megatonnes';
set(c, 'fontsize', 20); set(gca, 'fontsize', 20);
set(gcf, 'position', [0 0 1400 600]); 

%% Reachable structural capacity

% We remove the first entry of 'tmass' to obtain a vector where each entry
% represents the structural capacity in mass terms of the trap with
% correspodning index. (The removed element represents the region outside any
% spill region, so it has a trapping capacity of zero).
tmass = tmass(2:end);

cum_reachable = zeros(size(ta.traps));

% Adding the capacity of each trap to its own spill region and all downstream
% regions
for trap_ix = 1:max(ta.traps)

   region = ta.trap_regions==trap_ix;

   % Counting this trap's capacity towards all cells in its spill region
   cum_reachable(region) = cum_reachable(region) + tmass(trap_ix);

   visited_regions = trap_ix;

   % Computing contribution to cells associated with downstream traps
   downstream = find(ta.trap_adj(:,trap_ix));
   while ~isempty(downstream)

      region = ta.trap_regions == downstream(1);
      cum_reachable(region) = cum_reachable(region) + tmass(trap_ix);

      visited_regions = [visited_regions;downstream(1)];

      % downstream(1) has now been processed, so we remove it from the
      % downstream vector.
      downstream = [downstream(2:end); find(ta.trap_adj(:,downstream(1)))];
      downstream = setdiff(downstream, visited_regions);
   end
end

% we divide 'cum_reachable' by 1e9 to have the result in megatons
cum_reachable = cum_reachable / 1e9;

% Plot result (normal view)
figure
plotCellData(Gt, cum_reachable, 'edgealpha', 0.1);
view(90, 90); axis tight; colormap cool; 
c = colorbar; c.Label.String = 'Megatonnes';
set(c, 'fontsize', 20); set(gca, 'fontsize', 20);
set(gcf, 'position', [0 0 1400 600]); 

% Inclined view
view(80, 50);


%% Johansen
G = makeJohansenVEgrid();
Gt = topSurfaceGrid(G);

ta = trapAnalysis(Gt, false);

% Plot traps, rivers and spill regions (same code as for Utsira)

% construct color field identifying traps
trapfield = ta.traps;
trapfield(trapfield==0) = NaN;

% construct color field identifying rivers
river_field = NaN(size(ta.traps));
for i = 1:numel(ta.cell_lines)
   rivers = ta.cell_lines{i};
   for r = rivers
      river_field(r{:}) = 1;
   end
end

% Top view

% plot traps in solid color
clf; plotCellData(Gt, trapfield, 'edgecolor', 'none'); 

% plot regions in semi-transparent color
plotCellData(Gt, ta.trap_regions, 'facealpha', 0.5, 'edgecolor', 'none');

% add rivers
plotCellData(Gt, river_field, 'edgealpha', 0);

view(-43, 70); axis tight; colormap lines;
set(gca, 'fontsize', 20);
set(gcf, 'position', [0 0 1400 600]); camlight headlight; material dull
