%% Water-alternating Gas Injection in layer 10 of SPE10
% In this example, we simulate water-alternating gas injection (WAG), where
% we inject small volumes of water and solvent gas into the reservoir
% several cycles. The solvent gas mixes with the reservoir hydrocarbons
% according to amodified vesrion of the Todd-Longstaff mixing model, which
% treats the solvent gas as fourth pseudocomponent.
mrstModule add ad-blackoil ad-props spe10 mrst-gui solvent

df = get(0, 'defaultfigureposition');
close all

%% Set up grid and rock properties
% We pick up layer 10 of SPE10, and extract the grid and rock properties.
% We will define our own fluid based on a three-phase fluid with water, oil
% and gas.
fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 3, 0.4]*centi*poise, ...
                           'c'     , [1e-7, 1e-6, 1e-5]/barsa);
layers = 10;
[~, model4Ph, ~] = setupSPE10_AD('layers', layers);
G    = model4Ph.G;
rock = model4Ph.rock;

%% Define solvent properties         
% The solvent model we will use, treats solvent gas as a fourth phase,
% which is either miscible or immiscible with the oil and gas, depending on
% the fraction of solvent concentration to total gas concentration,
% $S_s/(S_g + S_s)$, and the pressure.

sOres_i = 0.38; % Immiscible residual oil saturation
sOres_m = 0.21; % Miscible residual oil saturation
fluid   = addSolventProperties(fluid, 'n'    , 2                   , ...
                                    'rho'    , 100*kilogram/meter^3, ...
                                    'mixPar' , 2/3                 , ...
                                    'mu'     , 0.5*centi*poise     , ...
                                    'sOres_i', sOres_i             , ...
                                    'sOres_m', sOres_m             , ...
                                    'c'      , 1e-5/barsa          );

%% Inspect fluid relperms
% In regions with only oil, reservoir gas and water, we have traditional
% black-oil behavior, whereas regions with only water, reservoir oil and
% solvent gas, the oil is completely miscible with the solvent. In this
% case, solvent gas mixes with formation oil, effectively altering the
% viscosities, densities and relperms according the the Todd-Longstaff
% model [1]. In the intermediate region, we interpolate between the two
% extrema depending on the solvent to total gas saturation fraction, and
% the pressure. We look at terneary plots for the hydrocarbon relperms of
% the hydrocardon phases when no water is present. Notice how the residual
% oil saturation (white line) reduces from $S_{or,i}$ (immiscible) to
% $S_{or,m}$ (miscible) with increasing solvent saturation.

figure('Position', [df(1:2), 1000, 400]);

n = 100;
[sO,sS] = meshgrid(linspace(0,1,n)', linspace(0,1,n)');
sW = 0; sG = 1-(sW+sO+sS);

ss  = linspace(0,1-sOres_m,100)';
b   = sOres_i + 2*ss - 1;
sg  = (-b + sqrt(b.^2 - 4*(ss.*(sOres_m + ss - 1))))/2;
sOr = 1 - sg - ss;

sW = sW(:); sO = sO(:); sG = sG(:); sS = sS(:);
[sWres, sOres , sSGres ]  = computeResidualSaturations(fluid, 0 , sG , sS );
[krW, krO, krG, krS]      = computeRelPermSolvent(fluid, 0, sW, sO, sG, sS, sWres, sOres, sSGres, 1);
sO = reshape(sO, n,n); sG = reshape(sG, n,n); sS = reshape(sS, n,n);

phName = {'O', 'G', 'S'};

for phNo = 1:3

    subplot(1,3,phNo)
    
    relpermName = ['kr', phName{phNo}];
    kr = reshape(eval(relpermName),[n,n]);
    kr(sG<0) = nan;

    [mapx, mapy] = ternaryAxis('names', {'S_g', 'S_s', 'S_o'});
    contourf(mapx(sG, sS, sO), mapy(sG,sS,sO), kr, 20, 'linecolor', 0.5.*[1,1,1])
    [mapx, mapy] = ternaryAxis('names', {'S_g', 'S_s', 'S_o'});
    plot(mapx(sg, ss, sOr), mapy(sg,ss,sOr), 'color', 0.99*[1,1,1], 'linewidth', 2)
    
    axis([0,1,0,sqrt(1-0.5^2)]); axis equal
    title(relpermName, 'position', [0.5,-0.2]); 
    
end
                                
%% Set up wells and injection schedule
% We consider the hello-world of reservoir simulation: A quarter-five spot
% pattern. The injector in the middle is operater at a fixed injection
% rate, while the producers in each corner is operated at a fixed
% bottom-hole pressure of 50 bars. During the first 2 years, we inject only
% water, while during the next two, we perform four cycles of equal periods
% of solvent gas and water injection. This is a common strategy to reduce
% the growth of viscous fingers associated with gas injection, and uphold a
% more favourable injected to reservoir fluid mobility ratio.

time = 1*year;                             % Time of the water injection
                                           % and WAG injection
rate = 0.15*sum(poreVolume(G, rock))/time;  % We inject a total of 0.2 PVI 
bhp = 50*barsa;                            % Producer bottom-hole pressure

producers = {1                                                     , ...
             G.cartDims(1)                                         , ...
             G.cartDims(1)*G.cartDims(2) - G.cartDims(1) + 1       , ...
             G.cartDims(1)*G.cartDims(2)                           };
injectors = {round(G.cartDims(1)*G.cartDims(2)/2 + G.cartDims(1)/2)};

W = [];
for pNo = 1:numel(producers)
    W = addWell(W, G, rock, producers{pNo}   , ...
                'type'  , 'bhp'              , ...
                'val'   , bhp                , ...
                'comp_i', [1,0,0,0]          , ...
                'sign'  , -1                 , ...
                'name'  , ['P', num2str(pNo)]);
end
for iNo = 1:numel(injectors)
    W = addWell(W, G, rock, injectors{iNo}     , ...
                'type'  , 'rate'               , ...
                'val'   , rate/numel(injectors), ...
                'comp_i', [1,0,0,0]            , ...
                'sign'  , 1                    , ...
                'name'  , ['I', num2str(iNo)]  );
end

dt          = 60*day; % Timestep size for water injection
dtWAG       = 45*day; % Timestep size for WAG cycles
useRampUp   = true;   % We use a few more steps each time we change the 
                      % well control to ease the nonlinear solution process
nCycles     = 4;      % Four WAG cycles
scheduleWAG = makeWAGschedule(W, nCycles, 'time'     , 2*time     , ...
                                          'dt'       , dtWAG    , ...
                                          'useRampup', useRampUp);

tvec                  = rampupTimesteps(time, dt);
control               = 2*ones(numel(tvec),1);
scheduleWAG.step.val     = [tvec; scheduleWAG.step.val; tvec];
scheduleWAG.step.control = [control; scheduleWAG.step.control; control];

%% Define model and initial state, and simulate results
% The reservoir is initially filled with a mixture of oil, reservoir gas
% and water. We set up a four-phase solvent model and simulate the
% schedule.
model4Ph = FourPhaseSolventModel(G, rock, fluid);
model4Ph.extraStateOutput = true;
sO = 0.5; sG = 0.4;
state0 = initResSol(G, 100*barsa, [1-sO-sG, sO, sG, 0]);
state0.wellSol = initWellSolAD(W, model4Ph, state0);

fn = getPlotAfterStepSolvent(state0, model4Ph, schedule, ...
                     'plotWell', true, 'plotReservoir', true);

[wellSolsWAG, statesWAG, reportsWAG] ...
    = simulateScheduleAD(state0, model4Ph, scheduleWAG, 'afterStepFn', fn);

%% Interactive visualization of results
% We use an interactive plotting tool to visualize reservoir results and
% the well solutions.
figure('Position', [df(1:2), 600,800])
plotToolbar(G, statesWAG); axis equal tight

%% Compare with standard water injection
% We compare the results with standard water injection to see the effects
% of using WAG.
for wNo = 1:numel(W)
    W(wNo).compi = [1,0,0];
end

tvec     = rampupTimesteps(4*time, dt);
schedule = simpleSchedule(tvec, 'W', W);
model3Ph = ThreePhaseBlackOilModel(G, rock, fluid);
model3Ph.extraStateOutput = true;
state0 = initResSol(G, 100*barsa, [1-sO-sG, sO, sG]);
state0.wellSol = initWellSolAD(W, model3Ph, state0);

[wellSol, states, reports] = simulateScheduleAD(state0, model3Ph, schedule);

%% Visulize well production curves
% We look at the well production curves to see the difference of the two
% injection strategies.
close all
tWAG = cumsum(scheduleWAG.step.val);
t    = cumsum(schedule.step.val);
plotWellSols({wellSolsWAG, wellSol},{tWAG, t})
