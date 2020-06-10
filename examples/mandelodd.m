clear all
close all


%% Include necessary modules
mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech mrst-gui
gravity off; % we do not consider gravity in this example

%% define grid
physdim = [40, 20] * meter;
resolution = [40, 20];
G = cartGrid(resolution, physdim);
G = computeGeometry(G);

%% Setup rock parameters

% flow parameters
perm = 100 * milli * darcy;
poro = 0.25;

% elastic parameters
young = 1 * giga * Pascal;
poisson = 0.3;
alpha = 1; % biot's coefficient

rock.poro  = poro  * ones(G.cells.num, 1);
rock.perm  = perm  * ones(G.cells.num, 1);
rock.alpha = alpha * ones(G.cells.num, 1);

top_force = 100 * mega * Pascal;

% boundary conditions will be set in the 'mech' object further down
mech.E = young * ones(G.cells.num, 1);
mech.nu = poisson * ones(G.cells.num, 1);
mech.load = @(x) 0*x;  % no body force applied

%% setup fluid parameters

muW = 0.89 * milli * Pascal / second;
rhoW = 1000 * kilogram / meter^3;
cW = 5.1e-10 * Pascal^-1;
pRef = 0;%1 * atm;

fluid = initSimpleADIFluid('phases' , 'W'  , ... % only water phase
                           'mu'     , muW  , ...
                           'rho'    , rhoW , ...
                           'c'      , cW   , ...
                           'pRef'   , pRef);

%% Set boundary conditions

% flow
bc = fluxside([], G, 'Ymin', 0); % no flow through bottom
bc = fluxside(bc, G, 'Ymax', 0); % no flow through top

bc = pside(bc, G, 'Xmin', pRef); % flow through left side
bc = pside(bc, G, 'Xmax', pRef); % flow through right side

% define lambda function to identify nodes for a set of faces
facenodes = ...
    @(f) unique(G.faces.nodes(mcolon(G.faces.nodePos(f), ...
                                     G.faces.nodePos(f+1)-1)));
% no displacement at bottom 
% @@ FIX THIS - lateral displacement should be allowed!

bfaces = getfield(pside([], G, 'Ymin', 1), 'face'); % bottom faces
bnodes = facenodes(bfaces); % bottom nodes
nbn = numel(bnodes); % number of bottom nodes

% no lateral displacement at sides
sfaces = [getfield(pside([], G, 'Xmin', 1), 'face'); ...
          getfield(pside([], G, 'Xmax', 1), 'face')]; % side faces
snodes = facenodes(sfaces); % side nodes
nsn = numel(snodes); % number of side nodes

disp_bc = struct('nodes', [bnodes; snodes], ...
                 'uu', zeros(nbn + nsn, 2), ...
                 'mask', [ones(nbn, 2); repmat([1, 0], nsn, 1)]);

% constant downward force on top
tfaces = getfield(pside([], G, 'Ymax', 1), 'face');

force_bc = struct('faces', tfaces, ...
                  'force', repmat([0, -1], numel(tfaces), 1) * top_force);

% combine displacement and force boundary conditions
mech.el_bc = struct('disp_bc', disp_bc, 'force_bc', force_bc);

%% setup model, schedule and initial state

% model
model = MechWaterModel(G, rock, fluid, mech);

model.FacilityModel = FacilityModel(model); % @@ hack - fix this!

% initial state
mech_unknowns = ~model.mechModel.operators.isdirdofs;

initState.pressure = pRef * ones(G.cells.num, 1);
initState.xd = zeros(nnz(mech_unknowns), 1);
initState = addDerivedQuantities(model.mechModel, initState);

% schedule
tsteps = 100;
duration = 10 * second;
schedule.step.val = duration/tsteps * ones(tsteps, 1);
schedule.step.control = ones(tsteps, 1); % no wells, so no controls
schedule.control = struct('W', [], 'bc', bc); % no wells

%% Run simulation
[~, states] = simulateScheduleAD(initState, model, schedule);


%% plotting
figure 

plotToolbar(G, states);

figure
ind = floor(resolution(1)/2);
pmid = cellfun(@(state) state.pressure(ind), states);
tt = cumsum(schedule.step.val);
plot(tt, pmid, '*');