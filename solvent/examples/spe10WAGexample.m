%% Water-alternating Gas Injection in layer 10 of SPE10
% In this example, we simulate water-alternating gas injection (WAG) in a
% layer of SPE10 model 2, where we inject small volumes of water and
% solvent gas into the reservoir several cycles. The solvent gas mixes with
% the reservoir hydrocarbons according to a modified version of the
% Todd-Longstaff mixing model, which treats the solvent gas as fourth
% pseudocomponent.
mrstModule add ad-blackoil ad-props spe10 mrst-gui solvent

df = get(0, 'defaultfigureposition');
close all

%% Set up grid and rock properties
% We pick up layer 10 of SPE10, and extract the grid and rock properties.
layers        = 10;
[~, model, ~] = setupSPE10_AD('layers', layers);
G             = model.G;
rock          = model.rock;

%% Set up fluid and add solvent gas properties
% We will define our own fluid based on a three-phase fluid with water, oil
% and gas. The solvent model we will use, treats solvent gas as a fourth
% pseudocomponent, which is either miscible or immiscible with the oil and
% gas, depending on the fraction of solvent concentration to total gas
% concentration, $S_s/(S_g + S_s)$, and the pressure. The model assumes one
% immiscible and one miscible residual saturation for the hydrocarbon
% phases, with $S_{\alpha r,i} > S_{\alpha r,m}$. In this example, we
% assume that the critical (residual) gas saturation is zero in both the
% immiscible and miscible case. The model uses a mixing paramter $\omega$
% that defines the degree of mixing, where (no mixing) $ = 0 <= \omega <= 1
% = $ (full mixing).
fluid = initSimpleADIFluid('n'     , [2, 2, 2]                        , ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG'                            , ...
                           'mu'    , [1, 3, 0.4]*centi*poise          , ...
                           'c'     , [1e-7, 1e-6, 1e-5]/barsa         );

sOr_i = 0.38; % Immiscible residual oil saturation
sOr_m = 0.21; % Miscible residual oil saturation

% We scale the relperms to the immiscible endpoints.
fluid.krW = coreyPhaseRelpermAD(2, 0, fluid.krG(1-sOr_i), sOr_i);
fluid.krG = coreyPhaseRelpermAD(2, 0, fluid.krG(1-sOr_i), sOr_i);     
[fluid.krO, fluid.krOW, fluid.krOG] = deal(coreyPhaseRelpermAD(2, sOr_i, fluid.krO(1), sOr_i));

fluid = addSolventProperties(fluid, 'rhoSS' , 100*kilogram/meter^3, ...
                                      'mixPar', 2/3                 , ...
                                      'muS'   , 0.5*centi*poise     , ...
                                      'sOr_m' , sOr_m               , ...
                                      'c'     , 1e-5/barsa          );

% We use dynamic end-point scaling to account for the changes in residual
% saturaitons. This is disabeled by default since it can be somewhat
% unstable.
model4Ph = BlackOilSolventModel(G, rock, fluid, 'dynamicEndPointScaling', true);

%% Inspect fluid relperms
% In regions with only oil, reservoir gas and water, we have traditional
% black-oil behavior, whereas regions with only water, reservoir oil and
% solvent gas, the oil is completely miscible with the solvent. In this
% case, solvent gas mixes with formation oil, effectively altering the
% viscosities, densities and relperms according the the Todd-Longstaff
% model. In the intermediate region, we interpolate between the two extrema
% depending on the solvent to total gas saturation fraction, and the
% pressure. We look at terneary plots for the hydrocarbon relperms of the
% hydrocardon phases when no water is present. Notice how the residual oil
% saturation (white line) reduces from $S_{or,i}$ (immiscible) to
% $S_{or,m}$ (miscible) with increasing solvent saturation.
plotSolventFluidProps(model4Ph, {'kr'}, {'O', 'G', 'S'});
                                
%% Set up wells and injection schedule
% We consider the hello-world of reservoir simulation: A quarter-five spot
% pattern. The injector in the middle is operater at a fixed injection
% rate, while the producers in each corner is operated at a fixed
% bottom-hole pressure of 50 bars. During the first half year, we inject
% only water, while during the next year, we perform four cycles of equal
% periods of solvent gas and water injection. Finally, we perform half a
% year of water injection. This is a common strategy to reduce the growth
% of viscous fingers associated with gas injection, and uphold a more
% favourable injected to reservoir fluid mobility ratio.

timeW   = 0.5*year;                           % Duration of water injection periods
timeWAG = 2*timeW;                            % Duration of WAG period
rate    = 0.1*sum(poreVolume(G, rock))/timeW; % We inject a total of 0.4 PVI 
bhp     = 50*barsa;                           % Producer bottom-hole pressure

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

dt1         = 60*day; % Timestep size for first water injection period
dt2         = 30*day; % Timestep size for WAG cycles and second water injection period
useRampUp   = true;   % We use a few more steps each time we change the 
                      % well control to ease the nonlinear solution process
nCycles     = 4;      % Four WAG cycles
scheduleWAG = makeWAGschedule(W, nCycles, 'time'     , timeWAG  , ...
                                          'dt'       , dt2      , ...
                                          'useRampup', useRampUp);

tvec1                    = rampupTimesteps(timeW, dt1);
tvec2                    = rampupTimesteps(timeW, dt2, 0);
scheduleWAG.step.val     = [tvec1; scheduleWAG.step.val; tvec2];
scheduleWAG.step.control = [2*ones(numel(tvec1),1)  ; ...
                            scheduleWAG.step.control; ...
                            2*ones(numel(tvec2),1) ];

%% Define model and initial state, and simulate results
% The reservoir is initially filled with a mixture of oil, reservoir gas
% and water. We set up a four-phase solvent model and simulate the
% schedule.
model4Ph.extraStateOutput = true;
sO = 0.5; sG = 0.4;
state0 = initResSol(G, 100*barsa, [1-sO-sG, sO, sG, 0]);
state0.wellSol = initWellSolAD(W, model4Ph, state0);

fn = getPlotAfterStepSolvent(state0, model4Ph, scheduleWAG, ...
                     'plotWell', true, 'plotReservoir', true);
% Use compiled, iterative linear solver with CPR preconditioning
mrstModule add linearsolvers
lsol = selectLinearSolverAD(model4Ph);
                 
[wellSolsWAG, statesWAG, reportsWAG] ...
    = simulateScheduleAD(state0, model4Ph, scheduleWAG, 'LinearSolver', lsol, 'afterStepFn', fn);

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

tvec     = rampupTimesteps(4*timeW, dt1);
schedule = simpleSchedule(tvec, 'W', W);
model3Ph = ThreePhaseBlackOilModel(G, rock, fluid);
model3Ph.extraStateOutput = true;
state0 = initResSol(G, 100*barsa, [1-sO-sG, sO, sG]);
state0.wellSol = initWellSolAD(W, model3Ph, state0);

lsol = selectLinearSolverAD(model3Ph);

[wellSol, states, reports] = simulateScheduleAD(state0, model3Ph, schedule, 'LinearSolver', lsol);

%% Visulize well production curves
% We look at the well production curves to see the difference of the two
% injection strategies.
close all
tWAG = cumsum(scheduleWAG.step.val);
t    = cumsum(schedule.step.val);
for sNo = 1:numel(wellSol)
    [wellSol{sNo}.qSs] = deal(0);
end
plotWellSols({wellSolsWAG, wellSol},{tWAG, t}, ...
             'datasetnames', {'WAG', 'Water flooding'});

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
