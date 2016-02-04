
opt.modelname                   = 'Utsirafm';
opt.coarse_level                = 5;
opt.max_num_wells               = 50; % updated to be a limit @@

opt.maximise_boundary_distance      = false; % if true, the following is not respected
opt.well_buffer_dist                = 1 * kilo * meter; % dist from edge of internal catchment
opt.well_buffer_dist_domain         = 5 * kilo * meter; % dist from edge of domain    
opt.well_buffer_dist_catchment      = 3 * kilo * meter; % dist from edge of external catchment
opt.pick_highest_pt                 = true; % Otherwise, furthest downslope is selected.
                                             % Both options are within constraints.
                                             

info                            = getSeaInfo(opt.modelname, 760);
opt.rhoW                        = info.water_density;
opt.sw                          = info.res_sat_wat;
opt.ref_temp                    = info.seafloor_temp + 273.15; % Kelvin
opt.ref_depth                   = info.seafloor_depth;
opt.temp_grad                   = info.temp_gradient;


%% Physical grid and rock
[Gt, rock2D, ~] = getFormationTopGrid(opt.modelname, opt.coarse_level);
if any(isnan(rock2D.perm))
    rock2D.perm = 500*milli*darcy * ones(Gt.cells.num,1);
end
if any(isnan(rock2D.poro))
    rock2D.poro = 0.25 * ones(Gt.cells.num,1); 
end

%% Spill-point analysis objec
ta = trapAnalysis(Gt, false);

%% CO2 property object
co2 = CO2props();
   
%% Testing well placement algorithm:
% [wc, qt] = pick_wellsites_test(Gt, rock2D, co2, ta, opt.max_num_wells, opt.rhoW, ...
%                         opt.sw, opt.ref_temp, opt.ref_depth, opt.temp_grad, ...
%                         opt.well_buffer_dist, opt.maximise_boundary_distance, ...
%                         opt.well_buffer_dist_domain, opt.pick_highest_pt, ...
%                         opt.well_buffer_dist_catchment, ...
%                         'inspectWellPlacement',true); 
                    
% [wc, qt] = pick_wellsites_array(Gt, rock2D, co2, ta, opt.rhoW, ...
%                         opt.sw, opt.ref_temp, opt.ref_depth, opt.temp_grad, ...
%                         1*kilo*meter, 1*kilo*meter, 2*kilo*meter, 2*kilo*meter, ...
%                         opt.max_num_wells, 'inspectWellPlacement',true);
                    
[wc, qt] = pick_wellsites_onePerTrapRegion(Gt, rock2D, co2, ta, opt.rhoW, ...
                        opt.sw, opt.ref_temp, opt.ref_depth, opt.temp_grad, ...
                        opt.max_num_wells, opt.pick_highest_pt, ...
                        opt.well_buffer_dist, opt.well_buffer_dist_domain, ...
                        opt.well_buffer_dist_catchment, ...
                        'inspectWellPlacement',true); 

%% Well placement inspection:
isteps = 50;
itime = 50 * year;
msteps = 10;
mtime = 10 * year;
opt.schedule = setSchedule(Gt, rock2D, wc, qt/info.rhoCref, isteps, ...
                                     itime, msteps , mtime, true, ...
                                     'minval', sqrt(eps));
cinx_inj = [opt.schedule.control(1).W.cells];

%figure(100); 
figure;
set(gcf,'Position',[2929 666 953 615])
clf

subplot(2,2,[1 3])
mapPlot(gcf, Gt, 'traps', ta.traps, ...
    'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
    'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
    'maplines', 20, 'wellcells', cinx_inj, 'well_numbering', true);
colorizeCatchmentRegions(Gt, ta);
plotGrid(Gt, 'FaceColor','none','EdgeAlpha',0.1)
axis equal tight off

subplot(2,2,2); bar(qt);
xlabel('well number', 'FontSize',16);
ylabel('mass to inject [kg]', 'FontSize',16); % i.e., mass capacity of empty traps along spill path

rates = [opt.schedule.control(1).W.val].*info.rhoCref.*(1*year)./10^9; % Mt/yr
subplot(2,2,4); bar(rates);
xlabel('well number', 'FontSize',16);
ylabel({'initial rate [Mt/yr]';['for ',num2str(convertTo(itime,year)),' yrs inj.']}, 'FontSize',16);


%sum(qt)/1e9/1e3