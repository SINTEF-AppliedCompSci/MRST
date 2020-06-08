
clear all
close all


%% load modules

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

%% Setup default options
opt = struct('cartDim'            , [100, 10], ...
             'L'                  , [100, 10], ...
             'fluid_model'        , 'oil water', ...
             'method'             , 'fully coupled', ...
             'nonlinearTolerance' , 1e-6, ...
             'verbose'            , false);

%% Setup grid


G = cartGrid(opt.cartDim, opt.L);
G = computeGeometry(G);

%% Setup rock parameters (for flow)

rock.perm = darcy*ones(G.cells.num, 1);
rock.poro = 0.3*ones(G.cells.num, 1);
% this is the Biot parameter which will be used to multiply the given dilatation
rock.alpha = ones(G.cells.num, 1);

%% Setup fluid parameters from SPE1

pRef = 270*barsa;
switch opt.fluid_model
  case 'blackoil'
    pth = getDatasetPath('spe1');
    fn  = fullfile(pth, 'BENCH_SPE1.DATA');
    deck = readEclipseDeck(fn);
    deck = convertDeckUnits(deck);
    fluid = initDeckADIFluid(deck);
    if isfield(fluid, 'pcOW')
        fluid = rmfield(fluid, 'pcOW');
    end
    if isfield(fluid, 'pcOG')
        fluid = rmfield(fluid, 'pcOG');
    end

    % Setup quadratic relative permeabilities, since SPE1 relperm are a bit rough.
    fluid.krW = @(s) s.^2;
    fluid.krG = @(s) s.^2;
    fluid.krOW = @(s) s.^2;
    fluid.krOG = @(s) s.^2;
    pRef = deck.PROPS.PVTW(1);

    modelargs = {'oil', true, 'gas', true, 'water', true};
    
  case {'oil water'}
    fluid = initSimpleADIFluid('phases', 'WO', 'mu', [1, 10]*centi*poise, ...
                               'n',  [1, 1], 'rho', [1000, 700]*kilogram/ ...
                               meter^3, 'c', 1e-10*[1, 1], 'cR', 4e-10, ...
                               'pRef', pRef);

    modelargs = {'oil', true, 'gas', false, 'water', true};
    
  case {'water'}
    fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                               1000*kilogram/meter^3, 'c', 1e-10, 'cR', ...
                               4e-10, 'pRef', pRef);
    
    modelargs = {'oil', false, 'gas', false, 'water', true};
  
  otherwise
    error('fluid_model  not recognized.');
end




%% Gravity
% The gravity in this option affects only the fluid behavior
gravity off;

fluidmodel = BiotFluidModel(G, rock, fluid, modelargs{:});


%% Setup wells
W = [];
refdepth = min(G.cells.centroids(:, G.griddim)) ;

ind = ceil(G.cartDims/2);
injcell = sub2ind(G.cartDims, ind(1), ind(2));

warning('off', 'RefDepth:BelowTopConnection');

W = addWell(W, G, rock, injcell, ...
            'Type'    , 'rate', ...
            'Val'     , 1e2/day, ...
            'Sign'    , 1,  ...
            'Comp_i'  , [0, 0, 1], ... % inject gas
            'Name'    , 'inj',  ...
            'refDepth', refdepth);

% W = addWell(W, G, rock, prodcell, ...
%             'Type'    ,'bhp', ...
%             'Val'     , pRef, ...
%             'Sign'    , -1,  ...
%             'Comp_i'  , [0, 1, 0], ... % one-phase test case
%             'Name'    , 'prod',  ...
%             'refDepth', refdepth);

switch opt.fluid_model
  case 'blackoil'
    W(1).compi = [0, 0, 1];
  case 'oil water'
    W(1).compi = [1 0];
    W(1).val   = 1e4/day;
  case 'water'
    W(1).compi = [1];
    W(1).val  = 1e-3/day;
  otherwise
    error('fluid_model not recognized.')
end

facilityModel = FacilityModel(fluidmodel);
facilityModel = facilityModel.setupWells(W);
model.FacilityModel = facilityModel;



%% Setup schedule
schedule.step.val     = [1*day*ones(1, 1); 5*day*ones(20, 1)];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control      = struct('W', W);

%% Setup initial state
clear initState;
initState.pressure = pRef*ones(G.cells.num, 1);
switch opt.fluid_model
  case 'blackoil'
    init_sat = [0, 1, 0];
    initState.rs  = 0.5*fluid.rsSat(initState.pressure);
  case 'oil water'
    init_sat = [0, 1];
  case 'water'
    init_sat = [1];
  otherwise
    error('fluid_model not recognized.')
end
nc = G.cells.num;
initState.s = ones(nc, 1)*init_sat;
initState.dilatation = 0.1*ones(nc, 1);

solver = NonLinearSolver('maxIterations', 100);
[wsol, states] = simulateScheduleAD(initState, fluidmodel, schedule, ...
                                    'nonlinearsolver', solver);

%% plot results
plotToolbar(G, states);
