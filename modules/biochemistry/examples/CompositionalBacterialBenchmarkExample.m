clear all;
mrstModule add compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on

%% Read the Eclipse deck file containing the simulation data
% Change input fil by UHS_BENCHMARK_RS_SALT.DATA for SALT EFFECTS
%deck = readEclipseDeck('/home/elyes/Documents/mrst-2023b/spe11-utils/deck_H2/UHS_benchmark/UHSDATA/UHS_BENCHMARK_RS.DATA');
deck = readEclipseDeck('/home/elyes/Documents/mrst-2023b/spe11-utils/TUC_UHS_Benchmark/Simulation Cases/UHS_Benchmark_HighH2.DATA');

%% Prepare simulation parameters and initial state
%[~, options, state0, model, schedule,compFluid, ~] = modified_uhs_benchmark_compositional(deck, 'bacteriamodel',false);
bacteriamodel = true;
require ad-core ad-props ad-blackoil
deck.PROPS.EOS ='PR';
deck = convertDeckUnits(deck);
compFluid = TableCompositionalMixture({'water', 'Methane','Nitrogen','CarbonDioxide','Ethane','Propane','n-butane','Hydrogen'}, ...
{'H2O', 'C1', 'N2', 'CO2', 'C2','C3','NC4','H2'});
eos =SoreideWhitsonEquationOfStateModel([], compFluid, 'sw');
[~, model, schedule, nls] = initEclipseProblemAD(deck,'getSchedule',true,'getInitialState', false);
% eos.PropertyModel.volumeShift = model.EOSModel.PropertyModel.volumeShift;
model.EOSModel =eos;
model.water = false;
T0 = deck.PROPS.TEMPVD{1}(2);
P0 = 82*barsa();
s0 = [0.2 0.8]; % hack to inforce residual saturation 
z0 = deck.PROPS.ZMFVD{1}(2:end);
z0 = [0.9723 0.00015  0.00405 0.020210 0.00272 0.0004 0.0002 0.000];
%z0 (end) = 0.0001;
% z0 (end-1) = 0.0001;
G = model.G;

%% Extend build up phase
[~, options,state0b, modelb, scheduleb, ~] = modified_uhs_benchmark(deck);
% Copy fields from scheduleb to schedule
% for i = 1:length(scheduleb.control)
% % scheduleb.control(i).W.cells = schedule.control(1).W.cells;
% % scheduleb.control(i).W.r = schedule.control(2).W.r;
% % scheduleb.control(i).W.dir = schedule.control(2).W.dir;
% % scheduleb.control(i).W.rR = schedule.control(2).W.rR;
% % scheduleb.control(i).W.WI = schedule.control(2).W.WI;
% % scheduleb.control(i).W.dZ = schedule.control(2).W.dZ;
% % scheduleb.control(i).W.heel = schedule.control(2).W.heel;
% % scheduleb.control(i).W.refDepth = schedule.control(1).W.refDepth;
% % scheduleb.control(i).W.defaulted = schedule.control(1).W.defaulted;
% % scheduleb.control(i).W.lims = schedule.control(1).W.lims;
% % scheduleb.control(i).W.cp = schedule.control(1).W.cp;
% scheduleb.control(i).W.components = schedule.control(1).W.components;
% % scheduleb.control(i).W.vfp_index = schedule.control(1).W.vfp_index;
% % scheduleb.control(i).W.cstatus = schedule.control(1).W.cstatus;
% % scheduleb.control(i).W.cell_origin = schedule.control(1).W.cell_origin;
% end
% schedule = scheduleb;
%%%
[rhow,rhog]=deal(model.fluid.rhoWS,0.081688* kilogram/meter^3); %density kilogram/meter^3;
[viscow,viscog]=deal(model.fluid.muWr,0.0094234*centi*poise);%viscosity
[cfw,cfg]=deal(model.fluid.cW,8.1533e-3/barsa); %compressibility
 %initialisation fluides incompressibles, Brooks-Corey relperm krw=(Sw)^nw
fluid=initSimpleADIFluid('phases', 'OG', 'mu',[viscow,viscog],...
                         'rho',[rhow,rhog],'pRef',0*barsa(),...
                         'c',[cfw,cfg],'n',[2,2],'smin',[0.2,0.05]);
fluid.krG = model.fluid.krG;
fluid.krO = model.fluid.krW;
fluid.krPts.g = model.fluid.krPts.g;
fluid.krPts.o = model.fluid.krPts.w;
% fluid.rhoOS = model.fluid.rhoWS;
% fluid.pvMultR = model.fluid.pvMultR;
% fluid.muO = model.fluid.muW;
% fluid.bO = model.fluid.bW;
fluid.pcOG = model.fluid.pcWG;

model.fluid = fluid;
deck = model.inputdata;
gravity reset on;
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);
if bacteriamodel
    compFluid = model.EOSModel.CompositionalMixture;
    arg = {model.G, model.rock, model.fluid, compFluid,...
        false, diagonal_backend, 'eos',model.EOSModel, 'oil', true, 'gas', true,... % water-oil system
    	'bacteriamodel', true,'diffusioneffect',false,'liquidPhase', 'O',...
        'vaporPhase', 'G'}; % water=liquid, gas=vapor
    model = BiochemistryModel(arg{:});
%     model.EOSModel =eos;
    model.outputFluxes = false;
    nbact0 = 1.0e8;
%     model.water = true;
     state0 = initCompositionalStateBacteria(model, P0.*ones(G.cells.num,1) , T0, s0, z0, nbact0, eos);
%      SimpleInitializationEquil;
%      state0.nbact = nbact0.*ones(G.cells.num,1);
%      state0.components( model.G.cells.centroids(:,3)>1210,:) = state0.components( model.G.cells.centroids(:,3)>1210,:).*0+ [0.9723 0.00015  0.00405 0.020210 0.00272 0.0004 0.0002 0];
%      state0.pressure = state0.pressure.*0 +(1.05*model.G.cells.centroids(:,3)./max(model.G.cells.centroids(:,3))).*78.*barsa;
%      state0.s = state0.s.*0 +(1.05*model.G.cells.centroids(:,3)./max(model.G.cells.centroids(:,3))).*1;
%      state0 = initCompositionalStateBacteria(model, state0.pressure , T0, s0, state0.components, nbact0, model.EOSModel);
%      state0.s = [sum(state0.s(:,1:2),2), state0.s(:,3)];
%      model.water = false;
else
    model.EOSModel =eos;
    state0 = initCompositionalState(model, P0.*ones(G.cells.num,1) , T0, s0, z0);
%     state0.components( model.G.cells.centroids(:,3)>1210,:) = state0.components( model.G.cells.centroids(:,3)>1210,:).*0+ [0.9723 0.00015  0.00405 0.020210 0.00272 0.0004 0.0002 0];
%     state0 = initCompositionalState(model, P0.*ones(G.cells.num,1) , T0, state0.components, z0);
end
% 
% 
% %===Ajout d'un terme source====================
%% ADD boundary condition
f = boundaryFaces(G);

%% Select lateral faces
f1 = (G.faces.normals(f, 3) > 2.48);
faces = f(~f1);

%% Calculate the displacement from a reference point (1140 units)
dx = bsxfun(@minus, G.faces.centroids(faces, :), 1140);

%% Calculate pressure drop due to gravity using omega and displacement
rhoOS = 9.98e+02;
dp = rhoOS .* (dx * reshape(gravity, [], 1));

%% Define the pressure at the boundary (init pressure plus the gravity-induced pressure drop)
pressure = state0.pressure(1) + dp;

%% Set the boundary condition with the computed pressure and saturation
bc = addBC([], faces, 'pressure', pressure, 'sat', state0.s(1,:));
bc.components = repmat(state0.components(1,:), numel(bc.face), 1);

for i=1:length(schedule.control)
    schedule.control(i).W.compi=[0, 1];
%     schedule.control(i).bc = bc;
end

%===Resolution pression/transport========================
% deltaT = T/niter;
% schedule = simpleSchedule(repmat(deltaT,1,niter),'bc', bc,'src', src,'W',W);
% nls = NonLinearSolver('useRelaxation', true);
%% Set up the linear and nonlinear solvers
lsolve = selectLinearSolverAD(model);                          % Select the linear solver for the model
nls = NonLinearSolver();                                       % Create a nonlinear solver object
nls.LinearSolver = lsolve;                                     % Assign the linear solver to the nonlinear solver

name = 'UHS_BENCHMARK_COMPOSITIONAL_BACT_TRUE_HIGH_H2_MOD';
%% Pack the simulation problem with the initial state, model, and schedule
% model = OverallCompositionCompositionalModel(G, model.rock,model.fluid,model.EOSModel, 'water', false);
% model.EOSModel.minimumComposition =1.0e-8;
 problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);
%[ws, states, reports] = simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);
%% Run the simulation
simulatePackedProblem(problem, 'restartStep',1);
%% gGet reservoir and well states
[ws,states] = getPackedSimulatorOutput(problem);
%% Compare with and without bectrial effects
problemNoBact = problem;
problemNoBact.BaseName = "UHS_BENCHMARK_COMPOSITIONAL_BACT_FALSE_HIGH_H2_MOD";
problemNoBact.OutputHandlers.states.dataDirectory= "/home/elyes/Documents/Projects/MRST/core/output/UHS_BENCHMARK_COMPOSITIONAL_BACT_FALSE_HIGH_H2_MOD";
problemNoBact.OutputHandlers.wellSols.dataDirectory= "/home/elyes/Documents/Projects/MRST/core/output/UHS_BENCHMARK_COMPOSITIONAL_BACT_FALSE_HIGH_H2_MOD";
problemNoBact.OutputHandlers.reports.dataDirectory= "/home/elyes/Documents/Projects/MRST/core/output/UHS_BENCHMARK_COMPOSITIONAL_BACT_FALSE_HIGH_H2_MOD";
[wsNoBact,statesNoBact] = getPackedSimulatorOutput(problemNoBact);

namecp = model.EOSModel.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indCO2= find(strcmp(namecp,'CO2'));
nT = numel(states);
for i = 1:nT
    totMassH2(i)=sum(states{i}.FlowProps.ComponentTotalMass{indH2}.*model.operators.pv);
    totMassH2(i)=sum(states{i}.FlowProps.ComponentTotalMass{indH2}.*model.operators.pv);
    totMassCO2(i)=sum(states{i}.FlowProps.ComponentTotalMass{indH2}.*model.operators.pv);

    totMassH2NoBact(i)=sum(statesNoBact{i}.FlowProps.ComponentTotalMass{indH2}.*model.operators.pv);
    totMassH2NoBact(i)=sum(statesNoBact{i}.FlowProps.ComponentTotalMass{indH2}.*model.operators.pv);
    totMassCO2NoBact(i)=sum(statesNoBact{i}.FlowProps.ComponentTotalMass{indH2}.*model.operators.pv);
end

