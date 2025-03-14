clear all;
mrstModule add compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on

%% Read the Eclipse deck file containing the simulation data
% Change input fil by UHS_BENCHMARK_RS_SALT.DATA for SALT EFFECTS
deck = readEclipseDeck('C:\Users\elyesa\OneDrive - SINTEF\Documents\Projects\TUC_UHS_Benchmark\Simulation Cases/UHS_Benchmark_HighH2.DATA');
%% Prepare simulation parameters and initial state
bacteriamodel = true;
deck.PROPS.EOS ='PR';
deck = convertDeckUnits(deck);
compFluid = TableCompositionalMixture({'water', 'Methane','Nitrogen','CarbonDioxide','Ethane','Propane','n-butane','Hydrogen'}, ...
{'H2O', 'C1', 'N2', 'CO2', 'C2','C3','NC4','H2'});
[~, model, schedule, nls] = initEclipseProblemAD(deck,'getSchedule',true,'getInitialState', false);

%% We use SW EOS
eos =SoreideWhitsonEquationOfStateModel([], compFluid, 'sw');
model.EOSModel =eos;
model.water = false;
%% Update fluid model to consider oil as water phase
[rhow,rhog]=deal(model.fluid.rhoWS,0.081688* kilogram/meter^3); %density kilogram/meter^3;
[viscow,viscog]=deal(model.fluid.muWr,0.0094234*centi*poise);%viscosity
[cfw,cfg]=deal(model.fluid.cW,8.1533e-3/barsa); %compressibility
% fluid
fluid=initSimpleADIFluid('phases', 'OG', 'mu',[viscow,viscog],...
                         'rho',[rhow,rhog],'pRef',0*barsa(),...
                         'c',[cfw,cfg],'n',[2,2],'smin',[0.2,0.05]);
fluid.krG = model.fluid.krG;
fluid.krO = model.fluid.krW;
fluid.krPts.g = model.fluid.krPts.g;
fluid.krPts.o = model.fluid.krPts.w;
fluid.pcOG = model.fluid.pcWG;
model.fluid = fluid;

%% Here we mimick the initialization from the benchmark (method 2 in Eclipse, not implemented in MRST)
T0 = deck.PROPS.TEMPVD{1}(2);
P0 = 82*barsa();
s0 = [0.2 0.8];
% z0 = deck.PROPS.ZMFVD{1}(2:end);
G = model.G;
InitCompositionForBenchmark();
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);
if bacteriamodel
    compFluid = model.EOSModel.CompositionalMixture;
    arg = {model.G, model.rock, model.fluid, compFluid,...
        false, diagonal_backend, 'eos',model.EOSModel, 'oil', true, 'gas', true,... % water-oil system
    	'bacteriamodel', bacteriamodel,'diffusioneffect',false,'liquidPhase', 'O',...
        'vaporPhase', 'G'}; % water=liquid, gas=vapor
    model = BiochemistryModel(arg{:});
    model.outputFluxes = false;
    nbact0 = 100;
    state0 = initCompositionalStateBacteria(model, P0.*ones(G.cells.num,1) , T0, s0, z0, nbact0, eos);
else
    model.EOSModel =eos;
    state0 = initCompositionalState(model, P0.*ones(G.cells.num,1) , T0, s0, z0);
end

for i=1:length(schedule.control)
    schedule.control(i).W.compi=[0, 1];
end

%% Set up the linear and nonlinear solvers
lsolve = selectLinearSolverAD(model);                          % Select the linear solver for the model
nls = NonLinearSolver();                                       % Create a nonlinear solver object
nls.LinearSolver = lsolve;                                     % Assign the linear solver to the nonlinear solver

name = 'UHS_BENCHMARK_COMPOSITIONAL_BACT_TRUE_HIGH_H2_MOD';
%% Pack the simulation problem with the initial state, model, and schedule 
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);
%% Run the simulation
simulatePackedProblem(problem, 'restartStep',1);
%% Compare with and without bectrial effects
modelNoBact = model;
modelNoBact.bacteriamodel = false;
state0NoBact = state0;
state0NoBact.nbact = 0;
nameNoBact = 'UHS_BENCHMARK_COMPOSITIONAL_BACT_FALSE_HIGH_H2_MOD';
problemNoBact = packSimulationProblem(state0, modelNoBact, schedule, nameNoBact, 'NonLinearSolver', nls);
%% Run the simulation
simulatePackedProblem(problemNoBact,  'restartStep',1);
%% Get reservoir and well states
[ws,states] = getPackedSimulatorOutput(problem);
[wsNoBact,statesNoBact] = getPackedSimulatorOutput(problemNoBact);

namecp = model.EOSModel.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indCO2= find(strcmp(namecp,'CO2'));
indCH4= find(strcmp(namecp,'C1'));
nT = numel(states);
% Initialize arrays to store total H2 mass
totalH2_bact = zeros(numel(states), 1);
totalH2_noBact = zeros(numel(statesNoBact), 1);

for i = 1:numel(states)
    % With bacterial effects
    totalH2_bact(i) = sum(states{i}.FlowProps.ComponentTotalMass{indH2});

    % Without bacterial effects
    totalH2_noBact(i) = sum(statesNoBact{i}.FlowProps.ComponentTotalMass{indH2});

    % With bacterial effects
    totalCO2_bact(i) = sum(states{i}.FlowProps.ComponentTotalMass{indCO2});

    % Without bacterial effects
    totalCO2_noBact(i) = sum(statesNoBact{i}.FlowProps.ComponentTotalMass{indCO2});
    % With bacterial effects
    totalCH4_bact(i) = sum(states{i}.FlowProps.ComponentTotalMass{indCH4});

    % Without bacterial effects
    totalCH4_noBact(i) = sum(statesNoBact{i}.FlowProps.ComponentTotalMass{indCH4});
end

%% Calculate percentage of H2 loss
H2_loss_percentage = ((totalH2_noBact - totalH2_bact) ./ totalH2_noBact) * 100;
%% Calculate percentage of CO2 loss
CO2_loss_percentage = ((totalCO2_noBact - totalCO2_bact) ./ totalCO2_noBact) * 100;
%% Calculate percentage of CH4 production
CH4_loss_percentage = ((totalCH4_bact - totalCH4_noBact) ./ totalCH4_noBact) * 100;

%% Display final H2 loss
fprintf('Total H2 loss due to bacterial effects: %.2f%%\n', H2_loss_percentage(end));