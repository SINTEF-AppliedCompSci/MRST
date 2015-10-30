%% explore Greater Snohvit Field
% located in the Hammerfest Basin Aquifer, comprised of Tubaen, Nordmela
% (low permeability), and Sto formations.

% The Sto formation is subject to trapping analysis and the injection
% scenario, thus we refine the dataset when creating the grdecl_st. Then,
% trapping analysis is done to reveal structural traps which correspond
% quite closely to the prospects and structural closures depicted in
% chapter 6 of the Compile CO2 storage Atlas, by NPD. In particular, the
% Greater Snohvit structural closure is extracted from the Sto formation,
% and used for an injection scenario.

% For now, the injection scenario occurs in the Sto formation part, however
% a model could be constructed which is made up of all three layers, and
% the injection occurs in the Tubaen (bottom).

% Injection originally was into the Tubaen formation, but in 2011, the well
% was perforated at a shallower level, such that CO2 was injected into the
% Sto formation (Hansen et al 2013, pg 3566).

% Approximate locations of wells:
% well 7121/4-1: (9.254e5, 7.988e6)
% well 7120/6-1: (9.182e5, 7.989e6)
% well 7121/4-2: (9.210e5, 7.990e6)
% inj well 7121/4F-2H: (9.225e5, 7.988e6)
%   - ED 50, UTM sone 34 (5.01998e5, 7.945754e6)
%   - obtained from http://factpages.npd.no/factpages/Default.aspx?culture=nb-no&nav1=field&nav2=PageView%7cProducing&nav3=2053062
mrstModule add libgeometry opm_gridprocessing
mrstModule add mrst-gui
mrstVerbose on
gravity on

% directory for saving figures
figDirName = 'HammerfestFigs_test';
mkdir(figDirName)


%% Set coarsening (N) or refining (R) value
%N = 1;
R = 1;


%% Get datasets and construct Grids (G, Gt), rock properties

[grdecls, rd, petroinfo] = getAtlasGrid({'Stofm','Nordmelafm','Tubaenfm'},'coarsening',3);

% Each grid is processed and the 2D (VE) rock properties are obtained.
for i = 1:numel(grdecls)

    grd = grdecls{i};

    G         = processGRDECL(grd);
    G         = mcomputeGeometry(G);
    rock      = grdecl2Rock(grd, G.cells.indexMap);
    rock.perm = convertFrom(rock.perm, milli*darcy);
    [Gt, G]   = topSurfaceGrid(G);

    % Even though the rock properties for Sto, Nordmela, and Tubaen are
    % already 2D (no discretization in vertical direction), we still use
    % averageRock to get rock2D which contains one column array of perm,
    % (i.e., the perm in the x-direction), as VE models do not handle
    % anisotropy.
    rock2D    = averageRock(rock, Gt);
    
    % Renaming of grids and rock
    if strcmpi(grd.name,'Stofm')
        [Gt_st, G_st, rock_st, rock2D_st] = deal(Gt, G, rock, rock2D);
    elseif strcmpi(grd.name,'Nordmelafm')
        [Gt_nd, G_nd, rock_nd, rock2D_nd] = deal(Gt, G, rock, rock2D);
    elseif strcmpi(grd.name,'Tubaenfm')
        [Gt_tu, G_tu, rock_tu, rock2D_tu] = deal(Gt, G, rock, rock2D);
    end
end



%% Visualization to assess grids and rock properties.

% This figure may be adjusted for desired look in poster/paper:
% -----------------------------------------------
figure; set(gcf,'Position',[1 1 1000 800])
% top surface grid of Sto
%plotGrid(Gt_st,'FaceColor','none','EdgeColor','y')
% three formation grids, colored differently
%plotGrid(G_st,'FaceColor',[1 .9 .9],'EdgeColor','none')
hp1 = plotGrid(G_tu,'FaceColor','c','EdgeColor','none');
hp2 = plotGrid(G_nd,'FaceColor','m','EdgeColor','none');
hp3 = plotGrid(G_st,'FaceColor',[1 .9 .9],'EdgeColor','none');
%plotCellData(Gt_st, Gt_st.cells.z,'EdgeColor','none')
light('Position',[-1 -1 1],'Style','infinite');lighting phong
view([-100 30]); grid; axis equal tight off
hl = legend([hp3,hp2,hp1],{'St\o','Nordmela',['Tub',char(229),'en']},'Location','West');
set(hl,'Fontsize',20)
set(gca,'DataAspect',[1 1 0.05]);
title('Hammerfest Basin Aquifer')
set(gca,'Fontsize',20)

export_fig(gcf,[figDirName '/' 'HammerfestBasin3D'], '-png','-transparent')
% -----------------------------------------------

% Get color bar limits of the data, for consistency in plotting of three
% formations. CAUTION: Inf (or NaN) should be ignored.
ntg_max = max([ max(rock2D_st.ntg(isfinite(rock2D_st.ntg))), ...
                max(rock2D_nd.ntg(isfinite(rock2D_nd.ntg))), ...
                max(rock2D_tu.ntg(isfinite(rock2D_tu.ntg)))  ]);
            
poro_max = max([ max(rock2D_st.poro(isfinite(rock2D_st.poro))), ...
                 max(rock2D_nd.poro(isfinite(rock2D_nd.poro))), ...
                 max(rock2D_tu.poro(isfinite(rock2D_tu.poro)))  ]);
            
perm_max = max([ max(rock2D_st.perm(isfinite(rock2D_st.perm))), ...
                 max(rock2D_nd.perm(isfinite(rock2D_nd.perm))), ...
                 max(rock2D_tu.perm(isfinite(rock2D_tu.perm)))  ]);
             
ntg_min = min([ min(rock2D_st.ntg(isfinite(rock2D_st.ntg))), ...
                min(rock2D_nd.ntg(isfinite(rock2D_nd.ntg))), ...
                min(rock2D_tu.ntg(isfinite(rock2D_tu.ntg)))  ]);
            
poro_min = min([ min(rock2D_st.poro(isfinite(rock2D_st.poro))), ...
                 min(rock2D_nd.poro(isfinite(rock2D_nd.poro))), ...
                 min(rock2D_tu.poro(isfinite(rock2D_tu.poro)))  ]);
            
perm_min = min([ min(rock2D_st.perm(isfinite(rock2D_st.perm))), ...
                 min(rock2D_nd.perm(isfinite(rock2D_nd.perm))), ...
                 min(rock2D_tu.perm(isfinite(rock2D_tu.perm)))  ]);
             
depth_max = max([ max(Gt_st.cells.z), max(Gt_nd.cells.z), max(Gt_tu.cells.z) ]);
depth_min = min([ min(Gt_st.cells.z), min(Gt_nd.cells.z), min(Gt_tu.cells.z) ]);

% The following figure can be compared with the plots shown in Compiled CO2
% Atlas, chp 6, pg 129.
figure; set(gcf,'Position',[1 1 1480 1011])

mySubPlot = @(myGt, myData) plotCellData(myGt, myData, 'EdgeColor','none');

% Sto formation
subplot(3,4,1)
mySubPlot(Gt_st, Gt_st.cells.z); colorbar; caxis([depth_min depth_max])
ylabel('St\o','Fontweight','bold')
title({'Top surface depth (meter)';' '})
subplot(3,4,2)
mySubPlot(Gt_st, rock2D_st.ntg); colorbar; caxis([ntg_min ntg_max]); axis off;
title({'net-to-gross';' '})
subplot(3,4,3)
mySubPlot(Gt_st, rock2D_st.poro); hcbp = colorbar; hcbp.Tag = 'por'; caxis([poro_min poro_max]); axis off;
title({'porosity';' '})
subplot(3,4,4)
mySubPlot(Gt_st, rock2D_st.perm./(milli*darcy)); colorbar; caxis([perm_min perm_max]./(milli*darcy)); axis off;
title({'permeability (mD)';' '})

% Nordmela
subplot(3,4,5)
mySubPlot(Gt_nd, Gt_nd.cells.z); colorbar; caxis([depth_min depth_max])
ylabel('Nordmela','Fontweight','bold')
subplot(3,4,6)
mySubPlot(Gt_nd, rock2D_nd.ntg); colorbar; caxis([ntg_min ntg_max]); axis off;
subplot(3,4,7)
plotCellData(Gt_nd, rock2D_nd.poro); hcbp = colorbar; hcbp.Tag = 'por'; caxis([poro_min poro_max]); axis off;
subplot(3,4,8)
plotCellData(Gt_nd, rock2D_nd.perm./(milli*darcy)); colorbar; caxis([perm_min perm_max]./(milli*darcy)); axis off;

% Tubaen
subplot(3,4,9)
plotCellData(Gt_tu, Gt_tu.cells.z); colorbar; caxis([depth_min depth_max])
ylabel(['Tub',char(229),'en'],'Fontweight','bold')
subplot(3,4,10)
plotCellData(Gt_tu, rock2D_tu.ntg); colorbar; caxis([ntg_min ntg_max]); axis off;
subplot(3,4,11)
mySubPlot(Gt_tu, rock2D_tu.poro); hcbp = colorbar; hcbp.Tag = 'por'; caxis([poro_min poro_max]); axis off;
subplot(3,4,12)
mySubPlot(Gt_tu, rock2D_tu.perm./(milli*darcy)); colorbar; caxis([perm_min perm_max]./(milli*darcy)); axis off;

% Note: some spots in formations contain which may have Inf or NaN values.

% Adjust axis, fonts
hfig = gcf;
set(findobj(hfig.Children,'Type','axes'),'Fontsize',16,'box','on')
axis(findobj(hfig.Children,'Type','axes'),'equal','tight')
set(findobj(hfig.Children,'Type','axes'),'XTick',[900000 940000 980000],'YTick',[7920000 7960000 8000000 8040000])

% adjust ticks in colorbars
hcbs = findobj(hfig.Children,'Type','colorbar');
for i = 1:numel(hcbs)
    currTicks = hcbs(i).Ticks;
    set(hcbs(i),'Ticks',[hcbs(i).Limits(1) currTicks hcbs(i).Limits(2)])
    
    %set(hcbs(i),'XTickLabel',num2str(get(hcbs(i),'XTick'),'%.1f'))
    
    cbt = get(hcbs(i), 'Ticks');
    if strcmpi(hcbs(i).Tag, 'por')
        cbtl = arrayfun(@(t) sprintf('%.2f', t), cbt, 'UniformOutput', false);
    else
        cbtl = arrayfun(@(t) sprintf('%.1f', t), cbt, 'UniformOutput', false);
    end
    set(hcbs(i), 'TickLabels', cbtl)

end
% save figure
export_fig(gcf, [figDirName '/' 'FmRockProperties'], '-png','-transparent')

%% Here, we replace any inf values with the field's mean value.
% however, a better approach is to use the surrounding values to ensure no
% spikes in data
rocks = [rock2D_st; rock2D_nd; rock2D_tu];
for i = 1:numel(rocks)
    rock2D = rocks(i);
    
    rock2D.poro(isinf(rock2D.poro)) = mean( rock2D.poro(~isinf(rock2D.poro)) );
    rock2D.perm(isinf(rock2D.perm)) = mean( rock2D.perm(~isinf(rock2D.perm)) );
    rock2D.ntg(isinf(rock2D.ntg))   = mean( rock2D.ntg(~isinf(rock2D.ntg)) );
    
    rocks(i) = rock2D;
end
rock2D_st = rocks(1);
rock2D_nd = rocks(2);
rock2D_tu = rocks(3);


%% Set up injection location

% Determine cell index of Snohvit injection location
% Physical coordinate (approx, see chp 6 pg 135 in Atlas):
var.wellXcoord      = 9.225e5;
var.wellYcoord      = 7.988e6;
dv = bsxfun(@minus, Gt_st.cells.centroids(:,1:2), [var.wellXcoord, var.wellYcoord]);
[v,i] = min(sum(dv.^2, 2));
wellCellIndex = i; % or Gt.cells.indexMap(i);
%[wellInd_i, wellInd_j] = ind2sub(Gt_st.cartDims, wellCellIndex);
var.wellCoord_x = Gt_st.cells.centroids(wellCellIndex,1);
var.wellCoord_y = Gt_st.cells.centroids(wellCellIndex,2);

% Note: well cell index corresponds to Sto formation, but should be the
% same for the other formations


%% Visualize side profile slices of grid to detect obvious faulted areas
sliceIndex = [wellCellIndex; wellCellIndex + 500; wellCellIndex + 1000];
ta_st       = trapAnalysis(Gt_st, false);
[ ~, ~ ]    = Grid_withSideProfiles( Gt_st, var.wellXcoord, var.wellYcoord, ...
    var.wellCoord_x, var.wellCoord_y, ta_st, 'sliceCellIndex', sliceIndex, ...
    'plotNorthwardSlices', true);


%% Do some trapping analysis
% We will consider the top layer of the Hammerfest Basin aquifer first,
% which is the Sto formation. Q: which trapping estimates do porosity and
% permeability impact? Where has an averaged value (for entire formation)
% been used in estimates? Now that porosity and permeability maps are
% available, how will they impact trapping estimations?

%interactiveTrapping(Gt_st, 'injpt', wellCellIndex)
interactiveTrapping('Stofm', 'injpt', wellCellIndex) % use formation name in order to load averaged perm and poro data using getAtlasGrid
view([-85 45]);
% then use icons in figure to toggle various plotting options.
% Note: interactiveTrapping sets the poro and perm to default values of 0.3
% and 300 if Gt is passed in. If the formation name is passed in, it will
% get poro and perm values (non-nan) from getAtlasGrid.

% re-size figure for clarity
set(gcf,'Position',[1 1 822 839])
view(2)
hfig = gcf;
set(findobj(hfig.Children,'Type','axes'),'FontSize',14)
set(findobj(hfig.Children,'Type','Legend'),'FontSize',14,'Orientation','vertical')
haxs = hfig.Children; % pie axes should be haxs(5)?
set(findobj(haxs(4).Children,'Type','Text'),'FontSize',14) % pie chart might be haxes(4)

% save figure
export_fig(gcf, [figDirName '/' 'interactiveTrapping2'], '-png','-transparent')


% run a VE simulation (do manually)

% adjust views/position of some subplots
hfig = gcf;
haxs = findobj(hfig.Children,'Type','axes');
haxs(3).View = [0  90]; % (height of co2-column)
haxs(1).View = [-92 8]; % (co2-sat, y-slice)
haxs(2).View = [-1  8]; % (co2-sat, x-slice)
haxs(5).Position = [0.05 0.5838 0.4942 0.3412];
haxs(4).Position = [0.5 0.5838 0.4 0.3412];

% adjust font size
set(findobj(hfig.Children,'Type','axes'),'FontSize',14)
set(findobj(hfig.Children,'Type','Legend'),'FontSize',14,'Position',[0.7606 0.7151 0.2029 0.2000])
set(findobj(haxs(4).Children,'Type','Text'),'FontSize',14)

% save figure
export_fig(gcf, [figDirName '/' 'interactiveTrapping_VEsim'], '-png','-transparent')



% Note that injection wells could be placed in other locations of the Sto
% formation, to exploit the full capacity of the structural traps. Also, it
% is quite clear by the outline of some structural traps where some of the
% NPD's prospects are location. Trapping analysis should find these
% prospect areas and display detail such as rock volume, pore volume,
% storage capacity, etc.

% could also use:
exploreCapacity( 'default_formation',  'Stofm',     ...
                 'grid_coarsening',     1,          ...
                 'seafloor_depth',      330*meter,  ...
                 'seafloor_temp',       4,          ...
                 'temp_gradient',       40           );
             
exploreSimulation(  'default_formation',  'Nordmelafm',     ...
                    'grid_coarsening',     3,          ...
                    'seafloor_depth',      330*meter,  ...
                    'seafloor_temp',       4,          ...
                    'temp_gradient',       40,         ...
                    'inj_time',            30 * year,  ...
                    'inj_steps',           30,         ...
                    'mig_time',            0 * year,   ...
                    'mig_steps',           0,          ...
                    'savefile',            'testsave'   );
                
% note: averaged perm and poro is used in exploreSimulation() since
% getFormationTopGrid() is called which returns rock2D data with averaged
% properties. 


%% Run injection scenario using entire Sto formation
[wellSols_st, states_st, sim_report_st, opt_st, var_st ] = runSnohvitInjectionScenario( Gt_st, rock2D_st );


%% Run LARGE-SCALE injection scenario (for Sto formation)
[wellSols_st_ls, states_st_ls, sim_report_st_ls, opt_st_ls, var_st_ls ] = ...
    runSnohvitInjectionScenario( Gt_st, rock2D_st, ...
    'inj_rate_MtperYr', 10, 'inj_time', 50*year, 'inj_steps', 50, ...
    'mig_time', 3000*year, 'mig_steps', 30);

% Analyze simulated injection/migration:
figure;
plotToolbar(G_st,states_st_ls)

% Top View plot of final year CO2-sat 
Years2plot          = convertTo(sim_report_st.ReservoirTime(end), year);
SimStartYear        = 0;
plume               = [];
plumeOutline_SatTol = (0.01/100); % adjust this value if patch error occurs
                                  % (which can happen when plotting a year
                                  % with minimal massCO2)
    
[ ~, ~ ] = subplotCO2simVsCO2obsData_withSideProfiles_basic( Years2plot, ...
    SimStartYear, plume, sim_report_st_ls, Gt_st, states_st_ls, var_st_ls.fluid, var_st_ls.model, ...
    [], [], var_st_ls.wellCoord_x, var_st_ls.wellCoord_y, ...
    ta_st, plumeOutline_SatTol);
title([num2str(Year2plot),' years after sim. start'])

%% Alternatively, run LARGE-SCALE injection scenario using exploreSimulation
% compare 1 vs 3 injection wells located in the main structural traps
exploreSimulation(  'default_formation',  'Stofm',     ...
                    'grid_coarsening',     1,          ...
                    'seafloor_depth',      330*meter,  ...
                    'seafloor_temp',       4,          ...
                    'temp_gradient',       40,         ...
                    'inj_time',            50 * year,  ...
                    'inj_steps',           10,         ...
                    'mig_time',            3000 * year,   ...
                    'mig_steps',           30           );


%% Run injection scenario using entire Tubaen formation
% then assess the pressure build-up which has been previously studied
[wellSols_tu, states_tu, sim_report_tu, opt_tu, var_tu ] = runSnohvitInjectionScenario( Gt_tu, rock2D_tu );


%% Do some analysis
bf = boundaryFaces(Gt_st);
plotRealVsDiscreteInjLoc(Gt_st, bf, opt_st.wellXcoord, opt_st.wellYcoord, var_st.wellCoord_x, var_st.wellCoord_y);
makeSideProfilePlots_CO2heights(G_st, Gt_st, var_st.wellCellIndex, states_st, var_st.fluid );

% can do other analysis of results. i.e.,
figure;
plotToolbar(G_st,states_st)

%% Plot CO2-sat for a specific year (and plot side profiles)
% If a long-time injection/migration was simulated, plot CO2-sat/profiles
Years2plot          = convertTo(sim_report_st.ReservoirTime(end), year);
SimStartYear        = 0;
plume               = [];
plumeOutline_SatTol = (0.01/100); % adjust this value if patch error occurs
                                  % (which can happen when plotting a year
                                  % with minimal massCO2)
    
[ ~, ~ ] = subplotCO2simVsCO2obsData_withSideProfiles_basic( Years2plot, ...
    SimStartYear, plume, sim_report_st, Gt_st, states_st, var_st.fluid, var_st.model, ...
    [], [], var_st.wellCoord_x, var_st.wellCoord_y, ...
    ta_st, plumeOutline_SatTol);

% save figure
export_fig(gcf, [figDirName '/' 'VEsim_Stofm'], '-png','-transparent')


% Add side profiles:
% specify XY coordinate(s) you want to draw EW-slice(s) through:
XYcoord1_ew = [var.wellXcoord, var.wellYcoord + 600];
XYcoords_ew = [XYcoord1_ew];
% get index of closest cell to your specified XY coordinate(s):
cellIndex = zeros(numel(XYcoords_ew(:,1)),1);
for p = 1:numel(XYcoords_ew(:,1))
    dv              = bsxfun(@minus, var_st.model.G.cells.centroids(:,1:2), XYcoords_ew(p,:));
    [v, cinx]       = min(sum(dv.^2, 2));   
    cellIndex(p,:)  = cinx;
end
inxPts_ew = [cellIndex];
% specify XY coordinate(s) you want to draw SN-slice(s) through:
XYcoord1_sn = [var.wellXcoord + 400, var.wellYcoord];
XYcoords_sn = [XYcoord1_sn];
% get index of closest cell to your specified XY coordinate(s):
cellIndex = zeros(numel(XYcoords_sn(:,1)),1);
for p = 1:numel(XYcoords_sn(:,1))
    dv              = bsxfun(@minus, var_st.model.G.cells.centroids(:,1:2), XYcoords_sn(p,:));
    [v, cinx]       = min(sum(dv.^2, 2));   
    cellIndex(p,:)  = cinx;
end
inxPts_sn = [cellIndex];
[ ~, ~ ] = subplotCO2simVsCO2obsData_withSideProfiles_basic( Years2plot, ...
    SimStartYear, plume, sim_report_st, Gt_st, states_st, var_st.fluid, var_st.model, ...
    [], [], var_st.wellCoord_x, var_st.wellCoord_y, ...
    ta_st, plumeOutline_SatTol, ...
    'sliceCellIndex_ew',inxPts_ew, 'sliceCellIndex_sn',inxPts_sn );

% save figure
export_fig(gcf, [figDirName '/' 'VEsim_Stofm_withSideProfiles'], '-png','-transparent')


% preliminary trapping anaylsis on G. Snohvit

%% Determine segregation speed
% To validate assumption of gravity segratation and use of VE modeling, the
% segregation speed should be faster than the time it takes for CO2 to
% reach the caprock from the injection depth. (CONFIRM).
% A conservative estimate involves calculating the smaller (slowest)
% segregation speed.
drho    = opt_st.water_density - opt_st.co2_density;
permv   = mean(rock2D_st.perm);
phi     = mean(rock2D_st.poro);
seg_spd = drho * permv * norm(gravity) / ( opt_st.co2_mu * phi ); % m/s
fprintf('\n\n The segregation speed is %d meters/second. \n\n', seg_spd)
r4validVE = opt_st.inj_rate / ( 2 * pi() * seg_spd * Gt_st.cells.H(var_st.wellCellIndex) ); % m
fprintf('\n\n Given the injection rate of %d m3/s, gravity segregation is a valid assumption in the region beyond a radius of %d meters from the injection point. \n\n', opt_st.inj_rate, r4validVE)
fprintf('\n\n Note that the grid resolution is %d meters. \n\n', (Gt_st.cells.centroids(1,1) - Gt_st.cells.centroids(2,1)) )

% % eqn 10 in Ghanbarnezhad et al, 2015 (but for pore space velocity, not
% % darcy velocity.
% krG = var.fluid.krG( max(states{end}.s(:,2)) );
% krW = var.fluid.krG( min(states{end}.s(:,1)) ); % use min for conservative estimate?
% seg_spd_incomp  = drho * permv * norm(gravity) / (( opt.co2_mu/krG + opt.water_mu/krW ) * phi); % m/s



%% Cut out the Greater Snohvit field (GSF) from Sto

% Trapping analysis, to get structures in Sto formation (one of which we
% will assume to be the Greater Snohvit field)
ta_st = trapAnalysis(Gt_st, false);

% Location of a point within the Greater Snohvit field
snohvit_x      = 9.225e5;
snohvit_y      = 7.988e6;
dv = bsxfun(@minus, Gt_st.cells.centroids(:,1:2), [snohvit_x, snohvit_y]);
[~,i] = min(sum(dv.^2, 2));
snohvit_cellIndex = i; % or Gt.cells.indexMap(i);

% IDs of Greater Snohvit:
SnohvitTrapID   = ta_st.traps(snohvit_cellIndex);

% Cut away rest of Sto grid except cells included in this structural trap
G_gsf       = removeCells(G_st, ta_st.traps ~= SnohvitTrapID);
G_gsf.name  = 'Greater Snohvit';

% Get top surface and use it to extract rock2D for GSF
[Gt_gsf, G_gsf] = topSurfaceGrid(G_gsf);

% Extract rock2D for GSF
%rock2D_gsf  = removeCells(rock2D_st, ta_st.traps ~= SnohvitTrapID);
for i = 1:Gt_gsf.cells.num
    rock2D_gsf.perm(i,1) = rock2D_st.perm( find( Gt_st.cells.indexMap == Gt_gsf.cells.indexMap(i) ) );
    rock2D_gsf.poro(i,1) = rock2D_st.poro( find( Gt_st.cells.indexMap == Gt_gsf.cells.indexMap(i) ) );
    rock2D_gsf.ntg(i,1)  = rock2D_st.ntg( find( Gt_st.cells.indexMap == Gt_gsf.cells.indexMap(i) ) );
end
% optional: visualize perm, poro, ntg data of GSF.

% Perform trapping analysis of Greater Snohvit field
ta_gsf = trapAnalysis(Gt_gsf, false);

% Visualize Greater Snohvit field
figure;

% subplot to show location of Greater Snohvit field within Sto formation
subplot(1,2,1)
plotGrid(G_st,'FaceColor', [1 .9 .9], 'EdgeAlpha', .05);
%plotFaces(Gt_st, bf, 'EdgeColor','k', 'FaceColor','none',  'LineWidth',3);
plotGrid(G_gsf, 'FaceColor','r', 'EdgeAlpha', .05)
axis equal tight off
light('Position',[-1 -1 1],'Style','infinite');lighting phong
set(gca,'DataAspect',[1 1 0.02]);
view([-100 55]);
legend('Sto formation','Greater Snohvit')

% subplot to focus on grid of Greater Snohvit
subplot(1,2,2)
plotGrid(G_gsf, 'FaceColor','r', 'EdgeAlpha', .05)
axis equal tight off
light('Position',[-1 -1 1],'Style','infinite');lighting phong
set(gca,'DataAspect',[1 1 0.02]);
%view([-100 55]);
view(2)

%
figure;
mapPlot(gcf, Gt_gsf, 'maplines',10)%'traps',ta_gsf.traps)
% approximate actual location
plot(wellXcoord,wellYcoord,'ok', ...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10)
% simulated injection location
plot(wellCoord_x,wellCoord_y,'xk', ...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10)

figure;
subplot(1,4,1)
plotGrid(G_gsf, 'FaceColor','none'); view(3);
subplot(1,4,2)
plotCellData(G_gsf, rock2D_gsf.ntg); colorbar; title('ntg'); axis equal tight
subplot(1,4,3)
plotCellData(G_gsf, rock2D_gsf.poro); colorbar; title('poro'); axis equal tight
subplot(1,4,4)
plotCellData(G_gsf, rock2D_gsf.perm./(milli*darcy)); colorbar; title('perm'); axis equal tight
        
%% Visualize side profile slices of grid to detect obvious faulted areas
% first need to get the well cell index in the Gt_gsf grid:
dv = bsxfun(@minus, Gt_gsf.cells.centroids(:,1:2), [wellXcoord, wellYcoord]);
[v,i] = min(sum(dv.^2, 2));
wellCellIndex = i;
wellCoord_x = Gt_gsf.cells.centroids(wellCellIndex,1);
wellCoord_y = Gt_gsf.cells.centroids(wellCellIndex,2);

% then call to plotting routine
sliceIndex = [wellCellIndex; wellCellIndex + 100; wellCellIndex + 200];
[ ~, ~ ]    = Grid_withSideProfiles( Gt_gsf, wellXcoord, wellYcoord, ...
    wellCoord_x, wellCoord_y, ta_gsf, 'sliceCellIndex', sliceIndex, ...
    'plotNorthwardSlices', true);

        
%% Perform injection scenario in Greater Snohvit field
% region belonging to Greater Snohvit Field: see NDP figures in Chapter 6
% of Compiled Atlas
[wellSols_gsf, states_gsf, sim_report_gsf, opt_gsf, var_gsf ] = runSnohvitInjectionScenario( Gt_gsf, rock2D_gsf );
       

%% Do some analysis
bf = boundaryFaces(Gt_gsf);
plotRealVsDiscreteInjLoc(Gt_gsf, bf, opt_gsf.wellXcoord, opt_gsf.wellYcoord, var_gsf.wellCoord_x, var_gsf.wellCoord_y);
makeSideProfilePlots_CO2heights(G_gsf, Gt_gsf, var_gsf.wellCellIndex, states_gsf, var_gsf.fluid );

% can do other analysis of results. i.e.,
figure;
plotToolbar(G_gsf,states_gsf)
% Note: it might be better to cut out a region larger than the GSF
% structural trap, since boundary conditions along the sub-domaexin may be
% impacting CO2's containment inside the trap. Allowing it to spill out of
% domain may be more realistic.


%% Do some storage capacity estimation, to compare with NPD (pg 136)
% Here, we tabulate the data from NPD and compare MRST-estimates for the
% Greater Snohvit field to it. "Snohvit2800m" is the pore volume in Sto,
% Nordmela, and Tubaen formations in the Snohvit and Snohvit North area
% down to 2800m (NPD, pg 134). "SnohvitCentralSto" is the pore volume in Sto
% only, surrounding the G segment where good communication exists.
% Notes:    netVol * avPoro = poreVol
%           rockVol * avNTG = netVol
%           storageCapacity (m3) / poreVol (m3) = storageEfficiency (%) -->
%           see Bachu 2015, eqn 1
%  (where   storageCapacity (kg) / rhoCO2 (kg/m3) = storageCapacity (m3) )

% Bachu et al 2014 gives the following: E (CO2 storage efficiency
% coefficient) = 2.7% for dolomite, 2% for limestone, 2.4% for sandstone.

NDP_data.Snohvit2800m = struct( 'storage_system',       'half open', ...
                                'rock_volume_m3',       89 * 1e9 * meter^3, ...
                                'net_volume_m3',        54 * 1e9 * meter^3, ...
                                'pore_volume_m3',       6.4 * 1e9 * meter^3, ...
                                'max_depth_m',          2800 * meter, ...
                                'min_depth_m',          2404 * meter, ...
                                'average_ntg',          0.6, ...
                                'average_poro',         0.12, ...
                                'average_perm',         150 * milli * darcy, ...
                                'storage_efficiency',   2, ...
                                'storage_capacity_kg',  90 * 1e9 * kilo * gram );
NDP_data.Snohvit2800m.name = 'Snohvit2800m';
                            
NDP_data.SnohvitCentralSto = struct( 'storage_system',  'half open', ...
                                'rock_volume_m3',       6.1 * 1e9 * meter^3, ...
                                'net_volume_m3',        4.8 * 1e9 * meter^3, ...
                                'pore_volume_m3',       0.680 * 1e9 * meter^3, ...
                                'max_depth_m',          2400 * meter, ...
                                'min_depth_m',          2320 * meter, ...
                                'average_ntg',          0.8, ...
                                'average_poro',         0.14, ...
                                'average_perm',         300 * milli * darcy, ...
                                'storage_efficiency',   5, ...
                                'storage_capacity_kg',  24 * 1e9 * kilo * gram );
NDP_data.SnohvitCentralSto.name = 'SnohvitCentralSto';
                            
rhoCO2 = 760 * kilogram / meter ^3;
MRST_est.GSF = struct( 'storage_system',  'tbd', ...
                                'rock_volume_m3',       sum( G_gsf.cells.volumes ) * meter^3, ...
                                'net_volume_m3',        sum( G_gsf.cells.volumes ) * meter^3 * mean(rock2D_gsf.ntg), ...
                                'pore_volume_m3',       sum( G_gsf.cells.volumes ) * meter^3 * mean(rock2D_gsf.ntg) * mean(rock2D_gsf.poro), ...
                                'max_depth_m',          max(Gt_gsf.cells.z) * meter, ...
                                'min_depth_m',          min(Gt_gsf.cells.z) * meter, ...
                                'average_ntg',          mean(rock2D_gsf.ntg), ...
                                'average_poro',         mean(rock2D_gsf.poro), ...
                                'average_perm',         mean(rock2D_gsf.perm) * meter^2 );
                            
test_rock2D_gsf         = rock2D_gsf;
test_rock2D_gsf.poro    = rock2D_gsf.poro.*rock2D_gsf.ntg;
[ capacityOutput ] = getTrappingCapacities(Gt_gsf, test_rock2D_gsf, ta_gsf, ...
        opt_gsf.co2_density, opt_gsf.water_density, opt_gsf.seafloor_temp, opt_gsf.seafloor_depth, ...
        opt_gsf.temp_gradient, 0, opt_gsf.res_sat_co2, opt_gsf.res_sat_wat, opt_gsf.dis_max);
    
MRST_est.GSF.storage_capacity_kg =  capacityOutput.tot_trap_capa_sum * 1e3 * 1e9 * kilo * gram;
MRST_est.GSF.storage_efficiency = (( MRST_est.GSF.storage_capacity_kg / rhoCO2 ) / MRST_est.GSF.pore_volume_m3 ) * 100;
MRST_est.GSF.name = 'GreaterSnohvitField';

% Show table of data comparison
res = {NDP_data.Snohvit2800m; NDP_data.SnohvitCentralSto; MRST_est.GSF};
fprintf('\n\nOverview of trapping capacity:\n')
fprintf('Note: Capacity may or maynot refer to structural trap.\n')
%fprintf('\n%-13s| System  | Rock Vol  | Net Vol | Pore Vol | Max depth | Min depth | Av NTG | Av Poro | Av Perm | \n', 'Name');


   fprintf('----------------------|-----------------------|----------------------|----------------------\n');
   fprintf('Name                  | %-22s|%-22s|%-19s|\n',        res{1}.name,                  res{2}.name,                res{3}.name);
   fprintf('----------------------|-----------------------|----------------------|----------------------\n');
   fprintf('Storage system        | %-22s|%-22s|%-19s|\n',        res{1}.storage_system,        res{2}.storage_system,      res{3}.storage_system);
   fprintf('Rock vol (m3)         | %4.2e              | %4.2e             | %4.2e          |\n',   res{1}.rock_volume_m3,        res{2}.rock_volume_m3,      res{3}.rock_volume_m3);
   fprintf('Net vol (m3)          | %4.2e              | %4.2e             | %4.2e          |\n',   res{1}.net_volume_m3,         res{2}.net_volume_m3,       res{3}.net_volume_m3);
   fprintf('Pore vol (m3)         | %4.2e              | %4.2e             | %4.2e          |\n',   res{1}.pore_volume_m3,        res{2}.pore_volume_m3,      res{3}.pore_volume_m3);
   fprintf('Max depth (m)         | %4.2e              | %4.2e             | %4.2e          |\n',   res{1}.max_depth_m,           res{2}.max_depth_m,         res{3}.max_depth_m);
   fprintf('Min depth (m)         | %4.2e              | %4.2e             | %4.2e          |\n',   res{1}.min_depth_m,           res{2}.min_depth_m,         res{3}.min_depth_m);
   fprintf('Av ntg                | %4.2e              | %4.2e             | %4.2e          |\n',   res{1}.average_ntg,           res{2}.average_ntg,         res{3}.average_ntg);
   fprintf('Av poro               | %4.2e              | %4.2e             | %4.2e          |\n',   res{1}.average_poro,          res{2}.average_poro,        res{3}.average_poro);
   fprintf('Av perm               | %4.2e              | %4.2e             | %4.2e          |\n',   res{1}.average_perm,          res{2}.average_perm,        res{3}.average_perm);
   fprintf('Storage efficiency    | %4.2e              | %4.2e             | %4.2e          |\n',   res{1}.storage_efficiency,    res{2}.storage_efficiency,  res{3}.storage_efficiency );
   fprintf('Storage capacity (kg) | %4.2e              | %4.2e             | %4.2e          |\n',   res{1}.storage_capacity_kg,   res{2}.storage_capacity_kg, res{3}.storage_capacity_kg );
   fprintf('----------------------|-----------------------|----------------------|----------------------\n');



