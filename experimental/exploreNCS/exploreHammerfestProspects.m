%% explore Prospects and Structural Closures in Hammerfest Basin Aquifer
% The NPD has tabulated storage capacity estimates for several regions that
% are found within the Hammerfest Basin aquifer, most notably in the Sto
% formation. These regions are called prospects or structural closures.
% Here, we extract the prospects and structural closures from the Sto
% formation, by assuming they correspond to the structural traps found by
% the MRST algorithm of trapAnalysis(). Two prospects are found with a
% structural closure, so to extract those prospects, an additional trapping
% analysis is performed on the structural closure.

% Note, prospect G and H are not quite captured by an entire structural
% trap, and may require other properties to extract them from the Sto
% formation grid (such as elevation or permeability properties).


mrstModule add libgeometry opm_gridprocessing

N = 1;

%% Get Sto dataset
% First get Sto formation dataset and rock properties, and then get top
% surface grid. Note: perm is in mD
grdecl_st = addPermPoroNtgData2grdecl( 'Stofm', 'refining',N );
G_st = processGRDECL(grdecl_st);
G_st = mcomputeGeometry(G_st);
rock_st = grdecl2Rock(grdecl_st, G_st.cells.indexMap);
rock_st.perm = convertFrom(rock_st.perm, milli*darcy);
[Gt_st, G_st] = topSurfaceGrid(G_st);
rock2D_st  = averageRock(rock_st, Gt_st);



%% Do Quick Trapping Analysis

% quick way to see the individual structural traps (note that resolution of
% grid will impact trap connectivity):
ta_st = trapAnalysis(Gt_st, false);
figure;
plotCellData(Gt_st,ta_st.traps,'EdgeColor','none')

% to see contours, and structural traps:
figure;
mapPlot(gcf,Gt_st,'traps',ta_st.traps)


%% Cut out the structural closures and prospects

[ Grids ] = getProspectTrapIDsAndGrids( G_st, Gt_st, ta_st );


% visualization of Sto with Prospects and Greater structural closures.
figure;
plotGrid(G_st,'FaceColor', [1 .9 .9], 'EdgeAlpha', .05);
%plotFaces(Gt_st, bf, 'EdgeColor','k', 'FaceColor','none',  'LineWidth',3);
plotGrid(G_gsf, 'FaceColor',[1 0 0], 'EdgeAlpha', .05)
plotGrid(G_albatross, 'FaceColor',[1 0.4 0], 'EdgeAlpha', .05)
plotGrid(G_askeladd, 'FaceColor',[1 0.4 0.4], 'EdgeAlpha', .05)

plotGrid(G_prospectC, 'FaceColor','b', 'EdgeAlpha', .05)
plotGrid(G_prospectD, 'FaceColor','g', 'EdgeAlpha', .05)
plotGrid(G_prospectE, 'FaceColor','y', 'EdgeAlpha', .05)
plotGrid(G_prospectF, 'FaceColor',[0.5 0.7 0.2], 'EdgeAlpha', .05)
plotGrid(G_prospectG, 'FaceColor','m', 'EdgeAlpha', .05)
plotGrid(G_prospectH, 'FaceColor','c', 'EdgeAlpha', .05)

light('Position',[-1 -1 1],'Style','infinite');lighting phong

axis equal tight off
%light('Position',[-1 -1 1],'Style','infinite');lighting phong
set(gca,'DataAspect',[1 1 0.02]);
view([-100 55]);
legend('Sto formation','Greater Snohvit','Greater Albatross','Greater Askeladd','Prospect C','Prospect D','Prospect E','Prospect F','Prospect G','Prospect H', 'Location','NorthEast')

% another visualization (base to make prospect comparisons)
figure;
plotCellData(Gt_st, ta_st.traps, 'EdgeColor','none','FaceColor',[1 0.9 0.9])
names = fieldnames(Grids);
for i = 1:numel(names)
    G = getfield(Grids, names{i});
    Gt = topSurfaceGrid(G);
    %plotCellData(Gt, 'EdgeColor','none') % could select various cell data
    %to plot, to make comparisons between prospects. (ultimately, want to
    %colorize by storage capacity!)
    plotFaces(Gt, boundaryFaces(Gt), 'Linewidth',3)
end
light('Position',[-1 -1 1],'Style','infinite');lighting phong



%% Get volume and capacity of prospects, for comparison to NPD (Atlas chp 6)

% A -- using extracted Grids
rhoCO2 = 700; % kg/m3
Gs = [G_gsf, G_albatross, G_askeladd, G_prospectC, G_prospectD, ...
             G_prospectE, G_prospectF, G_prospectG, G_prospectH];
res = cell(numel(Gs),1);
fprintf('------------------------------------------------\n');
for i=1:numel(Gs)
   %fprintf('Processing %s ... ', grdecls{i}.name);
   G       = Gs(i);
   [Gt, ~] = topSurfaceGrid(G);
   ta      = trapAnalysis(Gt, false);
   
   res{i}.name      = G.name;
   res{i}.cells     = Gt.cells.num;
   res{i}.zmin      = min(Gt.cells.z);
   res{i}.zmax      = max(Gt.cells.z);
   res{i}.volume    = sum(G.cells.volumes);
   res{i}.trapvols  = volumesOfTraps(Gt,ta);
   res{i}.capacity_m3 = sum(res{i}.trapvols);
   res{i}.capacity_Mt = res{i}.capacity_m3*rhoCO2/(1e9);
   fprintf('done\n');
end

% Show table of volumes
fprintf('\n\nOverview of (structural) trapping capacity (in cubic meters):\n')
fprintf('Note: Rock volume and structural capacity are bulk values.\n')
fprintf('\n%-20s| Grid Cells |  Min  |  Max  | Rock Volume | Capacity (m3) | Capacity (Mt) | Percent\n', 'Name');
fprintf('--------------------|------------|-------|-------|-------------|---------------|---------------|--------\n');
for i=1:numel(Gs)
   fprintf('%-20s|   %6d   | %4.0f  | %4.0f  |   %4.2e  |    %4.2e   |    %4.2e   | %5.2f \n',...
      res{i}.name, res{i}.cells, res{i}.zmin, res{i}.zmax, res{i}.volume, ...
      res{i}.capacity_m3, res{i}.capacity_Mt, res{i}.capacity_m3/res{i}.volume*100);
end
fprintf('--------------------|------------|-------|-------|-------------|---------------|---------------|--------\n');


% B -- using traps remaining in Sto formation
% could write function where StructTrapIDs are given, then calculations are
% performed on each of these structural traps. The calculations are:
%   - spill pt elevation, top elevation of trap
%   - structural traps rock volume, av ntg, net vol, av poro, pore vol, av
%   perm
%   - max structural trapping capacity (in m3 and Mt)
%   - adjusted structural trapping capacity (in m3 and Mt) using E (storage
%   efficiency)



[rho, mu, sr, sw]   = getValuesSPE134891(); % from Sleipner benchmark
rhoCref = rho(2);
opt = struct(   'seafloor_depth',   330 * meter, ...
                'seafloor_temp',    4, ...
                'temp_gradient',    40, ...
                'water_density',    rho(1), ...
                'co2_density',      rho(2), ...
                'water_mu',         mu(1), ...
                'co2_mu',           mu(2), ...
                'res_sat_wat',      sw, ...
                'res_sat_co2',      sr, ...
                'water_compr_val',  4.3e-5/barsa, ...
                'pvMult',           1e-5/barsa, ...
                'isDissOn',         false, ...
                'dis_max',          (53 * kilogram / meter^3) / rhoCref );

for i = 1:numel(names)
    
    G       = getfield(Grids, names{i});
    [Gt, ~] = topSurfaceGrid(G);
    trapID  = G.trapID;
    trapName = { regexprep(G.name,'[^\w'']','') };  % needs to be passed in as a cell
    

    capacities{i} = getTrappingCapacities_specificTraps( Gt_st, rock2D_st, ta_st, ...
        trapID, trapName, ...
        opt.co2_density, opt.water_density, opt.seafloor_temp, opt.seafloor_depth, ...
        opt.temp_gradient, 0, opt.res_sat_co2, opt.res_sat_wat, opt.dis_max);
    
    
end
    

Tcurrent        = capacities{1}.T;
ColOfAmpSym     = table( repmat('&',numel(Tcurrent.Properties.RowNames),1), 'RowNames',Tcurrent.Properties.RowNames, 'VariableNames',{'AmpSym'}); % to insert between each column for latex table format
ColOfSlashSym   = table( repmat('\\',numel(Tcurrent.Properties.RowNames),1), 'RowNames',Tcurrent.Properties.RowNames, 'VariableNames',{'SlashSym'}); % to insert in last column for latex table format
Tcurrent        = join(ColOfAmpSym, Tcurrent, 'Keys','RowNames');
for i = 2:numel(capacities)
    T2add = capacities{i}.T;
    T2add = join( ColOfAmpSym, T2add, 'Keys','RowNames');
    Tcurrent = join( Tcurrent, T2add, 'Keys','RowNames');
end
Tcurrent = join( Tcurrent, ColOfSlashSym, 'Keys','RowNames');


%% Colorize Prospects/Structural closures by storage capacity (Mt)
% another visualization (base to make prospect comparisons)
figure;
plotCellData(Gt_st, ta_st.traps, 'EdgeColor','none','FaceColor',[1 0.9 0.9])
names = fieldnames(Grids);
for i = 1:numel(names)
    G = getfield(Grids, names{i});
    Gt = topSurfaceGrid(G);
    
    % colorize by storage capacity
    % TODO: ensure Grid and Capacity is for same trapName
    plotCellData(Gt, repmat(capacities{i}.storageCap_Mt, [Gt.cells.num,1]), 'EdgeColor','none')
    
    % plot boundaries of regions
    plotFaces(Gt, boundaryFaces(Gt), 'Linewidth',3)
end
light('Position',[-1 -1 1],'Style','infinite');lighting phong
axis equal tight off
hcb = colorbar; set(hcb,'FontSize',14);
set(hcb.Label, 'String','Mt', 'FontSize',16, 'Rotation',0);
set(gca,'DataAspect',[1 1 0.02]);
view([-100 55]);


% using mapPlot to show contours --> TODO: prevent colors from being
% changed by mapPlot
figure;
plotCellData(Gt_st, ta_st.traps, 'EdgeColor','none','FaceColor',[1 0.9 0.9])
mapPlot(gcf, Gt_st);
names = fieldnames(Grids);
for i = 1:numel(names)
    G = getfield(Grids, names{i});
    Gt = topSurfaceGrid(G);
    
    % colorize by storage capacity
    % TODO: ensure Grid and Capacity is for same trapName
    plotCellData(Gt, repmat(capacities{i}.storageCap_Mt, [Gt.cells.num,1]), 'EdgeColor','none')
    
    % plot boundaries of regions
    plotFaces(Gt, boundaryFaces(Gt), 'Linewidth',3)
end
mapPlot(gcf, Gt_st);
%light('Position',[-1 -1 1],'Style','infinite');lighting phong


% plot to show all structrual traps identified by trapAnalysis, and then
% the boundaries of the prospect regions ontop to show comparison
figure
mapPlot(gcf, Gt_st, 'traps',ta_st.traps)
names = fieldnames(Grids);
for i = 1:numel(names)
    G = getfield(Grids, names{i});
    Gt = topSurfaceGrid(G);
    % plot boundaries of regions
    plotFaces(Gt, boundaryFaces(Gt), [1 0 0], 'Linewidth',3)
end
axis equal tight off


%% Use NPD's storage efficiency values to get reduced storage capacity estimates



% other things to try:
% - increase resolution of initial grid, to get traps again
% - compare between cell-based and node-based method of trapAnalysis
% - use lower CO2 density (700 kg/m3) as used by NPD calculations
% - use formula with constant CO2 density to compute storage capacity


