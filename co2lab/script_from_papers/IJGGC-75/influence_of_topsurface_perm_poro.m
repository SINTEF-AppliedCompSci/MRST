%% Influence of caprock elevation, permeability and porosity on
% plume migration and structural trapping capacity

% The following script generates some of the results and plots presented in
% Section 3 of the paper, more specifically:
%   * Figure 6 (plume outlines and trapping inventories given perturbations
%   in permeability, porosity, and caprock elevation)
%   * Table 1 (middle column plots showing evolution of Sorensen-Dice
%   coefficient)
%   * Table 2 (data in 2nd and 4th row corresponding to caprock
%   perturbation of +/-5 meters and porosity perturbation of +/-0.05)
%   * Figure 7 (Histograms of structural capacities given perturbations in
%   caprock elevation and porosity)

% Note: num_real is the number of realizations that will be generated for
% each of the uncertain parameters considered here (caprock, permeability,
% porosity). An injection/migration scenario is simulated for each model
% realization, thus the following script can take several minutes to an
% hour, depending on how many realizations are generated. Results presented
% in the paper were generated using 100 realizations, however we specify
% only 10 realizations here for speed.
num_real = 10;

moduleCheck('ad-core','co2lab','mrst-gui')

%% Base model case (no perturbations applied yet)

[Gt_base, rock_base, seainfo] = makeMyGeomodel();
fluid_base = makeMyFluid(Gt_base, rock_base, seainfo);
schedule_base = makeMySchedule(Gt_base, rock_base, fluid_base);
initState_base = makeMyInitialState(Gt_base, seainfo);
model_base = CO2VEBlackOilTypeModel(Gt_base, rock_base, fluid_base);


% Simulating injection and migration:
[wellSols_base, states_base] = simulateScheduleAD(initState_base, model_base, schedule_base);


% Estimate structural trapping capacity (static)
ta_base = trapAnalysis(Gt_base, true); % independent of rock
tot_strap_base = computeTotalStructuralTrapping(Gt_base, ta_base, rock_base, seainfo); % depends on rock
    

% Some post-processing
reports_base = makeReports(model_base.G, {initState_base, states_base{:}}, ...
                model_base.rock, model_base.fluid, schedule_base, ...
                [model_base.fluid.res_water, model_base.fluid.res_gas], ...
                ta_base, []);
            
Ma_base = forecastCurve(initState_base, states_base, model_base, ta_base);

base.Gt_base         = Gt_base;
base.ta_base         = ta_base;
base.rock_base       = rock_base;
base.seainfo         = seainfo;
base.schedule        = schedule_base;
base.reports_base    = reports_base;
base.Ma_base         = Ma_base;


% Plot trapping structure of caprock
figure
mapPlot(gcf,Gt_base,'traps',ta_base.traps,'rivers',ta_base.cell_lines);
colorizeCatchmentRegions(Gt_base,ta_base);
axis equal tight off


%% Perturbed caprock
% Perturbation levels of +/-1, +/-5, and +/-15 meters were used to generate
% plots shown in Table 1 of paper. Here, we consider +/-5 meters.

% Loop through "r" realizations of caprock model
rock = rock_base;
for r=1:num_real
    
    % Generate a perturbed model
    [Gt, ~] = perturbMyGeomodel(Gt_base, rock, ...
                    'perturbType','perturbCaprock', ...
                    'pert_interval',[-5 5]);
    Gt_all{r} = Gt;
    
    
    % Estimate structural trapping capacity (static)
    ta = trapAnalysis(Gt, true); % independent of rock
    tot_strap(r) = computeTotalStructuralTrapping(Gt, ta, rock, seainfo); % depends on rock
   
    
    % Simulation injection/migration scenario for each perturbed model
    fluid = makeMyFluid(Gt, rock, seainfo);
    schedule = makeMySchedule(Gt, rock, fluid);
    initState = makeMyInitialState(Gt, seainfo);
    model = CO2VEBlackOilTypeModel(Gt, rock, fluid);

    [wellSols, states] = simulateScheduleAD(initState, model, schedule);

    
    % Prepare output
    reports_all{r} = makeReports(Gt, {initState, states{:}}, rock, fluid, ...
                        schedule, [fluid.res_water, fluid.res_gas], ta, []);
    Ma_all{r} = forecastCurve(initState, states, model, ta);

end

perturbed.Gt_all      = Gt_all;
perturbed.rock_all    = rock;
perturbed.tot_strap   = tot_strap;
perturbed.Ma_all      = Ma_all;
perturbed.reports_all = reports_all;


% Analyze results:
% Here, we specify select years for plotting plume outlines and error bars
% in the trapping inventory.
postProcess( base, perturbed, 'year',[30, 400, 800, 1200, 1600, 2000, 2400, 2800, 3000] );


% Plot histogram of structural trapping capacities using 100 realizations
% (Note that plot shown in Figure 7 of paper was generated using 1000
% realizations)
for r=1:100
    
    % Generate a perturbed model
    [Gt, ~] = perturbMyGeomodel(Gt_base, rock, ...
                    'perturbType','perturbCaprock', ...
                    'pert_interval',[-5 5]);
    
    % Estimate structural trapping capacity (static)
    ta = trapAnalysis(Gt, true); % independent of rock
    tot_strap(r) = computeTotalStructuralTrapping(Gt, ta, rock, seainfo); % depends on rock

end
plotSimplePDFCDF(tot_strap, 'numBins',15)


%% Perturbed permeability
% Perturbation levels of +/-0.02, +/-0.05, and +/-0.10 were used to
% generate plots shown in Table 1 of paper. Here, we consider +/-0.05.
% (Note that permeability perturbations are created by first perturbing
% porosity and then using a porosity-permeability model. As such the
% pert_interval specified here refers to a porosity perturbation. See
% paper, Section 2.3, for more explanation)

% Loop through "r" realizations of permeability model
Gt = Gt_base;
ta = ta_base;
for r=1:num_real
    
    % Generate a perturbed model
    [~, rock] = perturbMyGeomodel(Gt, rock_base, ...
                    'perturbType','perturbPerm', ...
                    'pert_interval',[-0.05 0.05]);
    rock_all{r} = rock;
    
    
    % Note: Structural trapping capacity (static) is not
    % dependent on permeability, thus we don't estimate it here.
   
    
    % Simulation injection/migration scenario for each perturbed model
    fluid = makeMyFluid(Gt, rock, seainfo);
    schedule = makeMySchedule(Gt, rock, fluid);
    initState = makeMyInitialState(Gt, seainfo);
    model = CO2VEBlackOilTypeModel(Gt, rock, fluid);

    [wellSols, states] = simulateScheduleAD(initState, model, schedule);

    
    % Prepare output
    reports_all{r} = makeReports(Gt, {initState, states{:}}, rock, fluid, ...
                        schedule, [fluid.res_water, fluid.res_gas], ta, []);
    Ma_all{r} = forecastCurve(initState, states, model, ta);

end

perturbed.Gt_all      = Gt;
perturbed.rock_all    = rock_all;
perturbed.tot_strap   = [];
perturbed.Ma_all      = Ma_all;
perturbed.reports_all = reports_all;


% Analyze results:
% Here, we specify select years for plotting plume outlines and error bars
% in the trapping inventory.
postProcess( base, perturbed, 'year',[30, 400, 800, 1200, 1600, 2000, 2400, 2800, 3000] );


%% Perturbed porosity
% Perturbation levels of +/-0.02, +/-0.05, and +/-0.10 were used to
% generate plots shown in Table 1 of paper. Here, we consider +/-0.05.

% Loop through "r" realizations of porosity model
Gt = Gt_base;
ta = ta_base;
for r=1:num_real
    
    % Generate a perturbed model
    [~, rock] = perturbMyGeomodel(Gt, rock_base, ...
                    'perturbType','perturbPoro', ...
                    'pert_interval',[-0.05 0.05]);
    rock_all{r} = rock;
    
    
    % Estimate structural trapping capacity (static)
    tot_strap(r) = computeTotalStructuralTrapping(Gt, ta, rock, seainfo); % depends on rock
   
    
    % Simulation injection/migration scenario for each perturbed model
    fluid = makeMyFluid(Gt, rock, seainfo);
    schedule = makeMySchedule(Gt, rock, fluid);
    initState = makeMyInitialState(Gt, seainfo);
    model = CO2VEBlackOilTypeModel(Gt, rock, fluid);

    [wellSols, states] = simulateScheduleAD(initState, model, schedule);

    
    % Prepare output
    reports_all{r} = makeReports(Gt, {initState, states{:}}, rock, fluid, ...
                        schedule, [fluid.res_water, fluid.res_gas], ta, []);
    Ma_all{r} = forecastCurve(initState, states, model, ta);

end

perturbed.Gt_all      = Gt;
perturbed.rock_all    = rock_all;
perturbed.tot_strap   = tot_strap;
perturbed.Ma_all      = Ma_all;
perturbed.reports_all = reports_all;


% Analyze results:
% Here, we specify select years for plotting plume outlines and error bars
% in the trapping inventory.
postProcess( base, perturbed, 'year',[30, 400, 800, 1200, 1600, 2000, 2400, 2800, 3000] );


% Plot histogram of structural trapping capacities using 100 realizations
% (Note that plot shown in Figure 7 of paper was generated using 1000
% realizations)
for r=1:100
    
    % Generate a perturbed model
    [~, rock] = perturbMyGeomodel(Gt, rock_base, ...
                    'perturbType','perturbPoro', ...
                    'pert_interval',[-0.05 0.05]);
   
    % Estimate structural trapping capacity (static)
    tot_strap(r) = computeTotalStructuralTrapping(Gt, ta, rock, seainfo); % depends on rock
   
end
plotSimplePDFCDF(tot_strap, 'numBins',15)