%% Add necessary modules

mrstModule add ad-blackoil ad-props % AD-moduels
mrstModule add mrst-gui
mrstModule add spe10                % SPE10 dataset
mrstModule add solvent              % Solvent model

gravity reset on

df = get(0, 'defaultfigureposition');
close all

%% Set up grid and rock properties
% We pick up layer 10 of SPE10, and extract the grid and rock properties.
% We will define our own fluid based on a three-phasespe10 fluid with water, oil
% and gas
fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [2, 3, 0.4]*centi*poise, ...
                           'c'     , [1e-7, 1e-6, 1e-5]/barsa);
layers = 10;
[~, model, ~] = setupSPE10_AD('layers', layers);
G    = model.G;
rock = model.rock;

%% Define solvent properties         
% The solvent model we will use, treats solvent gas as a fourth phase,
% which is either miscible or immiscible with the oil and gas, depending on
% the fraction of solvent concentration to total gas concentration,
% $S_s/(S_g + S_s)$, and the pressure.

sOres_i = 0.36; % Immiscible residual oil saturation
sOres_m = 0.12; % Miscible residual oil saturation
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
% case, olvent gas mixes with formation oil, effectively altering the
% viscosities, densities and relperms according the the Todd-Longstaff
% model [1]. In the intermediate region, we interpolate between the two
% extrema depending on the solvent to total gas saturation fraction, and
% the pressure. We look at terneary plots for the hydrocarbon relperms of
% the hydrocardon phases when no water is present. Notice how the residual
% oil saturation (white line) from $S_{or,i}$ (immiscible) to $S_{or,m}$
% (miscible) with increasing solvent saturation.

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
% bottom-hole pressure of 50 bars. We will simulate a popular injection
% strategy called water-alternating gas (WAG) injection over a period of
% four years. During the first 2 years, we inject only water, while the
% next two, we perform four cycles of equal periods of solvent gas and
% water injection. This is a common strategy to reduce the growth of
% viscous fingers associated with gas injection, and uphold a more
% favourable injected to reservoir fluid mobility ratio.

time = 1*year;                             % Time of the water injection and
                                           % WAG injection.
rate = 0.15*sum(poreVolume(G, rock))/time; % We inject a total of 0.3 PVI 
bhp = 50*barsa;                            % Producer bhp

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

dtWAG       = 15*day; % Timestep size
useRampUp   = true;   % We use a few more steps each time we change the 
                      % well control to ease the nonlinear solution process
nCycles     = 4;      % Four WAG cycles
scheduleWAG = makeWAGschedule(W, nCycles, 'time'     , time     , ...
                                          'dt'       , dtWAG    , ...
                                          'useRampup', useRampUp);
schedule = scheduleWAG;

dt                    = 30*day;
tvec                  = rampupTimesteps(time, dt);
control               = 2*ones(numel(tvec),1);
schedule.step.val     = [tvec; scheduleWAG.step.val];
schedule.step.control = [control; scheduleWAG.step.control];

%% Define model and initial state, and simulate results
% The reservoir is initially filled with a mixture of oil, reservoir gas
% and water.
model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;
sO = 0.75; sG = 0.2;
state0 = initResSol(G, 100*barsa, [1-sO-sG, sO, sG, 0]);
state0.wellSol = initWellSolAD(W, model, state0);

fn = getPlotAfterStepSolvent(state0, model, schedule, ...
                     'plotWell', true, 'plotReservoir', true);

[ws, states, reports] = simulateScheduleAD(state0, model, schedule, 'afterStepFn', fn);

%%

plotToolbar(G, states)
% plotWellSols(ws)

%%

