moduleCheck co2lab

%% Load Statfjord grid and perform spill point analysis
coarsening_level = 1; % grid downsampling factor (1 = no downsampling)
Gt = getFormationTopGrid('Statfjordfm', coarsening_level);
ta = trapAnalysis(Gt, false);

%% Show traps

figure
plotCellData(Gt, ta.traps, 'edgealpha', 0.1);
view(0, 90); axis tight; colormap lines;
set(gcf, 'position', [10 10 500 800]);

figure
plotCellData(Gt, ta.traps, 'edgealpha', 0.1);
view(290, 60); axis tight; colormap lines;
set(gcf, 'position', [520 10 1200 800]);


%% Show spill regions
trapfield = ta.traps;
trapfield(trapfield==0) = NaN;

figure;
plotCellData(Gt, trapfield, 'edgecolor', 'none');
plotCellData(Gt, ta.trap_regions, 'facealpha', 0.3, 'edgecolor', 'none');
view(0, 90); axis tight; colormap lines;
set(gcf, 'position', [10 10 500 800]);

figure;
plotCellData(Gt, trapfield, 'edgecolor', 'none');
plotCellData(Gt, ta.trap_regions, 'facealpha', 0.3, 'edgealpha', 0.1);
view(290, 60); axis tight; colormap lines;
set(gcf, 'position', [520 10 1200 800]);

%% Show size distributions in terms of area

% spill regions
region_cellcount = ...
    diff(find(diff([sort(ta.trap_regions); max(ta.trap_regions)+1])));
bar(sort(region_cellcount, 'descend'));
hold on;

% traps
trap_cellcount = diff(find(diff([sort(ta.traps); max(ta.traps)+1])));
bar(sort(trap_cellcount, 'descend'), 'r');

%% Show the largest spill region with associated trap

% Find index of trap with the largest spill region
[~, ix] = max(region_cellcount);

% Visualize the trap along with the spill region
figure;
plotGrid(topSurfaceGrid(extractSubgrid(Gt.parent, find(ta.traps==ix))), ...
         'facecolor', 'r', ...
         'edgealpha', 0.1)
plotGrid(topSurfaceGrid(extractSubgrid(Gt.parent, find(ta.trap_regions==ix))), ...
         'facecolor','r', ...
         'facealpha', 0.1, ...
         'edgealpha', 0.1);
view(-64, 36);
set(gcf, 'position', [10 10 1000 700]);


%% Show rivers
river_field = NaN(size(ta.traps));
for i = 1:numel(ta.cell_lines)
   rivers = ta.cell_lines{i};
   for r = rivers
      river_field(r{:}) = 1;
   end
end

figure;
plotCellData(Gt, river_field, 'edgealpha', 0.1); 
plotCellData(Gt, trapfield, 'edgecolor', 'none');
view(0, 90); axis tight; colormap lines;
set(gcf, 'position', [10 10 500 800]);

figure;
plotCellData(Gt, river_field, 'edgealpha', 0.1); 
plotCellData(Gt, trapfield, 'edgecolor', 'none');
view(290, 60); axis tight; colormap lines;
set(gcf, 'position', [520 10 1200 800]);


%% Diagram / map plot
h = figure;
mapPlot(h, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines);
set(gcf, 'position', [10 10 500 800]);

%% Compute trap volumes
depths = [0; ta.trap_z];

h = min(depths(ta.traps+1), Gt.cells.H);

cell_tvols = h .* Gt.cells.volumes;

tvols = accumarray(ta.traps+1, cell_tvols);
assert(tvols(1)==0);
tvols = tvols(2:end);

bar(sort(tvols, 'descend'))

%% Compute and plot trap capacity in mass terms

% We use the following assumptions for the parameters
porosity = 0.1071; % Porosity value for Statfjord from NPD
seafloor_temp = 7 + 273.15; % degrees Kelvin
temp_grad = 35; % degrees increase per kilometer depth
rho_brine = 1000; % brine density (kilogram per cubic meter)

% Ensure that gravity is not zero
gravity on;

% computing temperature field
T = seafloor_temp + temp_grad .* Gt.cells.z/1000;

% computing hydrostatic pressure
P = 1 * atm + rho_brine * norm(gravity) * Gt.cells.z;

% Making stripped-down fluid object that only contains CO2 density function
fluid = addSampledFluidProperties(struct, 'G');

cell_tmass = cell_tvols .* porosity .* fluid.rhoG(P, T);

tmass = accumarray(ta.traps+1, cell_tmass);
assert(tmass(1)==0);

figure
plotCellData(Gt, tmass(ta.traps+1), 'edgealpha', 0.1);
view(0, 90); axis tight; colormap cool; colorbar;
set(gcf, 'position', [10 10 500 800]);

figure
plotCellData(Gt, tmass(ta.traps+1), 'edgealpha', 0.1);
view(290, 60); axis tight; colormap cool; colorbar;
set(gcf, 'position', [520 10 1200 800]);


%% Compute reachable structural capacity

tmass = tmass(2:end); % remove the region spilling out of the domain
cum_reachable = zeros(size(ta.traps));

% Adding the capacity of each trap to its own spill region and all downstream
% regions
for trap_ix = 1:max(ta.traps)
   
   region = ta.trap_regions==trap_ix;
   cum_reachable(region) = cum_reachable(region) + tmass(trap_ix);
   
   visited_regions = trap_ix;
   
   downstream = find(ta.trap_adj(:,trap_ix));
   while ~isempty(downstream)
      
      region = ta.trap_regions == downstream(1);
      cum_reachable(region) = cum_reachable(region) + tmass(trap_ix);
      
      visited_regions = [visited_regions;downstream(1)];
      
      downstream = [downstream(2:end); find(ta.trap_adj(:,downstream(1)))];
      downstream = setdiff(downstream, visited_regions);
   end
end

% we divide 'cum_reachable' by 1e12 to have the result in gigatons
cum_reachable = cum_reachable / 1e12;

figure
plotCellData(Gt, cum_reachable, 'edgealpha', 0.1);
view(0, 90); axis tight; colormap cool; colorbar;
set(gcf, 'position', [10 10 500 800]);

figure
plotCellData(Gt, cum_reachable, 'edgealpha', 0.1);
view(290, 60); axis tight; colormap cool; colorbar;
set(gcf, 'position', [520 10 1200 800]);

   
   
   
   