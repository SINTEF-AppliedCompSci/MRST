%% Influence of caprock elevation, permeability and porosity on
% plume migration and structural trapping capacity

% The following script generates some of the results and plots presented in
% Section 3 of the paper (see README.txt), more specifically:
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
% hour, depending on how many realizations are generated. Simulation
% results presented in the paper were generated using 100 realizations,
% however we specify the number of realizations here (num_real) to be 10 in
% order to reduce computing time. Also, structural capacity estimates given
% perturbations in caprock and porosity were generated using 1000
% realizations, however we specify only 100 here (num_real_strap).
num_real = 10;
num_real_strap = 100;

mrstModule add co2lab
moduleCheck('ad-core')

%% Base model case (no perturbations applied yet)
% The geomodel we consider is the southern portion of the Utsira aquifer,
% located in the North Sea. This dataset is available as part of the
% Norwegian Petroleum Directorate's CO2 Storage Atlas. More details of this
% aquifer can be found here:
%       http://www.npd.no/en/Publications/Reports/CO2-Storage-Atlas-/
% The publically-downloadable Utsira geomodel is comprised of top and
% bottom surfaces, and an average permeability and porosity is reported in
% the published Atlas. Since heterogeneous rock properties are not
% available, we create plausible fields using a perm-phi and phi-depth
% model that we derive from another aquifer's geomodel.

% Using a single well, we simulate the injection of ~720 Mt of CO2 over a
% period of 30 years, and then a migration period of close to 3000 years.

[Gt_base, rock_base, seainfo] = makeMyGeomodel();
fluid_base = makeMyFluid(Gt_base, rock_base, seainfo);
schedule_base = makeMySchedule(Gt_base, rock_base, fluid_base);
initState_base = makeMyInitialState(Gt_base, seainfo);
model_base = CO2VEBlackOilTypeModel(Gt_base, rock_base, fluid_base);


% Simulating injection and migration:
[wellSols_base, states_base] = simulateScheduleAD(initState_base, ...
    model_base, schedule_base);


% Estimate structural trapping capacity (static)
ta_base = trapAnalysis(Gt_base, true); % independent of rock
tot_strap_base = computeTotalStructuralTrapping(Gt_base, ta_base, ...
    rock_base, seainfo); % depends on rock
    

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
% Now, we consider the influence of the caprock elevations on structural
% trapping capacity and plume migration by perturbing elevations using
% Guassian fields. Perturbation levels of +/-1, +/-5, and +/-15 meters were
% used to generate results shown in Table 1 of the paper. Here, we consider
% a perturbation level of +/-5 meters.

% Using same rock properties as base model
rock = rock_base;

% Loop through "r" realizations of caprock model
pertCaprock.Gt_all = cell(1,num_real);
pertCaprock.reports_all = cell(1,num_real);
pertCaprock.Ma_all = cell(1,num_real);
for r=1:num_real
    
    % Generate a perturbed model
    [Gt, ~] = perturbMyGeomodel(Gt_base, rock, ...
                    'perturbType','perturbCaprock', ...
                    'pert_interval',[-5 5]);
    pertCaprock.Gt_all{1,r} = Gt;

    
    % Simulation injection/migration scenario for each perturbed model
    fluid = makeMyFluid(Gt, rock, seainfo);
    schedule = makeMySchedule(Gt, rock, fluid);
    initState = makeMyInitialState(Gt, seainfo);
    model = CO2VEBlackOilTypeModel(Gt, rock, fluid);

    [wellSols, states] = simulateScheduleAD(initState, model, schedule);

    
    % Trapping structure is dependent on caprock elevations
    ta = trapAnalysis(Gt, true);
    
    
    % Prepare output
    pertCaprock.reports_all{r} = makeReports(Gt, {initState, states{:}}, ...
        rock, fluid, schedule, [fluid.res_water, fluid.res_gas], ta, []);
    pertCaprock.Ma_all{r} = forecastCurve(initState, states, model, ta);

end
pertCaprock.rock_all = rock;


% Analyze results:
% Here, we specify certain years for plotting plume outlines and error bars
% in the trapping inventory.
postProcess( base, pertCaprock, 'year', ...
    [30, 400, 800, 1200, 1600, 2000, 2400, 2800, 3000] );


% Now we create 100 realizations of the geomodel and plot a histogram of
% the cooresponding structural trapping capacities. Note that plot shown in
% Figure 7 of paper was generated using 1000 realizations.
pertCaprock.tot_strap = zeros(1,num_real_strap);
for r=1:num_real_strap
    
    % Generate a perturbed model
    [Gt, ~] = perturbMyGeomodel(Gt_base, rock, ...
                    'perturbType','perturbCaprock', ...
                    'pert_interval',[-5 5]);
    
    % Estimate structural trapping capacity (static)
    ta = trapAnalysis(Gt, true); % dependent on elevations
    pertCaprock.tot_strap(r) = computeTotalStructuralTrapping(Gt, ta, ...
        rock, seainfo);

end
plotSimplePDFCDF(pertCaprock.tot_strap, 'numBins',15)


%% Perturbed permeability
% Next, we consider the influence of permeability. Perturbation levels of
% +/-0.02, +/-0.05, and +/-0.10 were used to generate plots shown in Table
% 1 of paper. These levels actually refer to porosity perturbations, since
% permeability perturbations are created by first perturbing porosity and
% then using a porosity-permeability model. (See paper, Section 2.3, for
% more explanation.) Here, we consider a perturbation level of +/-0.05.
% Notice that we do not calculate structural trapping capacity for a set of
% permeability realizations here; this is because trapping capacity is a
% static estimate and is not dependent on permeability.

% Grid and trapping structure is not dependent on porosity
Gt = Gt_base;
ta = ta_base;

% Loop through "r" realizations of permeability model
pertPerm.rock_all = cell(1,num_real);
pertPerm.reports_all = cell(1,num_real);
pertPerm.Ma_all = cell(1,num_real);
for r=1:num_real
    
    % Generate a perturbed model
    [~, rock] = perturbMyGeomodel(Gt, rock_base, ...
                    'perturbType','perturbPerm', ...
                    'pert_interval',[-0.05 0.05]);
    pertPerm.rock_all{r} = rock;
    
    
    % Simulation injection/migration scenario for each perturbed model
    fluid = makeMyFluid(Gt, rock, seainfo);
    schedule = makeMySchedule(Gt, rock, fluid);
    initState = makeMyInitialState(Gt, seainfo);
    model = CO2VEBlackOilTypeModel(Gt, rock, fluid);

    [wellSols, states] = simulateScheduleAD(initState, model, schedule);

    
    % Prepare output
    pertPerm.reports_all{r} = makeReports(Gt, {initState, states{:}}, ...
        rock, fluid, schedule, [fluid.res_water, fluid.res_gas], ta, []);
    pertPerm.Ma_all{r} = forecastCurve(initState, states, model, ta);

end
pertPerm.Gt_all = Gt;


% Analyze results:
% Here, we specify certain years for plotting plume outlines and error bars
% in the trapping inventory.
postProcess( base, pertPerm, 'year', ...
    [30, 400, 800, 1200, 1600, 2000, 2400, 2800, 3000] );


%% Perturbed porosity
% Finally, we consider the influence of porosity. Perturbation levels of
% +/-0.02, +/-0.05, and +/-0.10 were used to generate plots shown in Table
% 1 of paper. Here, we consider perturbation level of +/-0.05.

% Grid and trapping structure is not dependent on porosity
Gt = Gt_base;
ta = ta_base;

% Loop through "r" realizations of porosity model
pertPoro.rock_all = cell(1,num_real);
pertPoro.reports_all = cell(1,num_real);
pertPoro.Ma_all = cell(1,num_real);
for r=1:num_real
    
    % Generate a perturbed model
    [~, rock] = perturbMyGeomodel(Gt, rock_base, ...
                    'perturbType','perturbPoro', ...
                    'pert_interval',[-0.05 0.05]);
    pertPoro.rock_all{r} = rock;
    
    
    % Simulation injection/migration scenario for each perturbed model
    fluid = makeMyFluid(Gt, rock, seainfo);
    schedule = makeMySchedule(Gt, rock, fluid);
    initState = makeMyInitialState(Gt, seainfo);
    model = CO2VEBlackOilTypeModel(Gt, rock, fluid);

    [wellSols, states] = simulateScheduleAD(initState, model, schedule);

    
    % Prepare output
    pertPoro.reports_all{r} = makeReports(Gt, {initState, states{:}}, ...
        rock, fluid, schedule, [fluid.res_water, fluid.res_gas], ta, []);
    pertPoro.Ma_all{r} = forecastCurve(initState, states, model, ta);

end
pertPoro.Gt_all = Gt;


% Analyze results:
% Here, we specify certain years for plotting plume outlines and error bars
% in the trapping inventory.
postProcess( base, pertPoro, 'year', ...
    [30, 400, 800, 1200, 1600, 2000, 2400, 2800, 3000] );


% Now we create 100 realizations of the geomodel and plot a histogram of
% the cooresponding structural trapping capacities. Note that plot shown in
% Figure 7 of paper was generated using 1000 realizations.
pertPoro.tot_strap = zeros(1,num_real_strap);
for r=1:num_real_strap
    
    % Generate a perturbed model
    [~, rock] = perturbMyGeomodel(Gt, rock_base, ...
                    'perturbType','perturbPoro', ...
                    'pert_interval',[-0.05 0.05]);
   
    % Estimate structural trapping capacity (static)
    pertPoro.tot_strap(r) = computeTotalStructuralTrapping(Gt, ta, rock, ...
        seainfo);
   
end
plotSimplePDFCDF(pertPoro.tot_strap, 'numBins',15)