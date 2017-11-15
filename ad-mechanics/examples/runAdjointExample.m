%% Computation of Adjoints for Lift values
% In this example, we demonstrate how one can use an adjoint simulation to
% compute gradients (sensitivities).
%
% The implementation and computational strategy used to obtain adjoints is
% independent of the actual system considered and a similar approach can be
% used for other systems and setups with minor (and obvious) modifications.
%

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

%% Setup geometry

% 2D cartesian domain
% todo: try some irregular geometry

cartDim = [100, 10];
L       = [100, 10];
G = cartGrid(cartDim, L);
G = computeGeometry(G);


%% Setup fluid

opt.fluid_model = 'oil water';
pRef = 100*barsa;
switch opt.fluid_model
  case 'oil water'
    init_sat = [0, 1];
  otherwise
    error('fluid_model not recognized.')
end

fluid = initSimpleADIFluid('phases', 'WO', 'mu', [1, 10]*centi*poise, 'n', ...
                           [1, 1], 'rho', [1000, 700]*kilogram/ meter^2, 'c', ...
                           1e-10*[1, 1], 'cR', 4e-10, 'pRef', pRef);


%% Setup rock parameters (for flow)

rock.perm = darcy*ones(G.cells.num, 1);
rock.poro = 0.3*ones(G.cells.num, 1);


%% Setup material parameters for Biot and mechanics

% same as 2D case
E          = 1*giga*Pascal; % Young's module
nu         = 0.3;           % Poisson's ratio
alpha      = 1;             % Biot's coefficient
% Convert global properties to cell values
E          = repmat(E, G.cells.num, 1);
nu         = repmat(nu, G.cells.num, 1);
rock.alpha = repmat(alpha, G.cells.num, 1);


%% Setup boundary conditions for mechanics (no displacement)
% zero displacement at bottom, left and right sides. Given pressure at the top.

% Gather the bc faces (left, bottom and right)
dummyval = 100; % We use pside to recover face at bottom, we use a dummy
                % value for pressure in this function.
bc = pside([], G, 'Xmin', dummyval);
bc = pside(bc, G, 'Xmax', dummyval);
bc = pside(bc, G, 'Ymin', dummyval);
indfacebc = bc.face;

% Get the nodes that belong to the bc faces.
facetonode = accumarray([G.faces.nodes, rldecode((1 : G.faces.num)', ...
                                                 diff(G.faces.nodePos))], ...
                        ones(numel(G.faces.nodes), 1), [G.nodes.num, ...
                    G.faces.num]);
isbcface = zeros(G.faces.num, 1);
isbcface(indfacebc) = 1;
bcnodes  = find(facetonode*isbcface);
nn       = numel(bcnodes);
u        = zeros(nn, G.griddim);
m        = ones(nn,  G.griddim);
disp_bc  = struct('nodes', bcnodes, 'uu', u, 'mask', m);

% Set a given pressure at the top face
dummyval = 100; % We use pside to recover face at bottom, we use a dummy
                % value for pressure in this function.
bc = pside([], G, 'Ymax', dummyval);
sidefaces = bc.face;
signcoef = (G.faces.neighbors(sidefaces, 1) == 0) - (G.faces.neighbors(sidefaces, ...
                                                  2) == 0);
n = bsxfun(@times, G.faces.normals(sidefaces, :), signcoef./ ...
           G.faces.areas(sidefaces));
force = bsxfun(@times, n, pRef);
force_bc = struct('faces', sidefaces, 'force', force);

% Construct the bc structure for mechanics
el_bc = struct('disp_bc' , disp_bc, 'force_bc', force_bc);


%% Setup volumetric load for mechanics
% In this example we do not impose any volumetric force
loadfun = @(x) (0*x);


%% Gather all the mechanical parameters in a struct
mech = struct('E', E, 'nu', nu, 'el_bc', el_bc, 'load', loadfun);


%% Set gravity off
gravity off


%% Setup model

switch opt.fluid_model
  case 'blackoil'
    error('not yet implemented')
  case 'oil water'
    model = MechOilWaterModel(G, rock, fluid, mech, 'verbose', true);
  case 'water'
  otherwise
    error('fluid_model not recognized.')
end


%% Set up initial reservoir state
clear initState;
initState.pressure = pRef*ones(G.cells.num, 1);
switch opt.fluid_model
  case 'blackoil'
    error('not yet implemented')
    init_sat = [0, 1, 0];
    initState.rs  = 0.5*fluid.rsSat(initState.pressure);
  case 'oil water'
    init_sat = [0, 1];
  case 'water'
    error('not yet implemented')
    init_sat = [1];
  otherwise
    error('fluid_model not recognized.')
end
initState.s  = ones(G.cells.num, 1)*init_sat;
initState.xd = zeros(nnz(~model.mechModel.operators.isdirdofs), 1);
initState    = computeInitDisp(model, initState, [], 'pressure', pRef* ...
                               ones(G.cells.num, 1));
initState    = addDerivedQuantities(model.mechModel, initState);


%% Setup the wells

nx = G.cartDims(1);
ny = G.cartDims(2);
W = [];
% injection wells
wellopt = {'type', 'rate', 'val', 1*meter^3/day, 'Sign', 1, 'comp_i', [1, 0]};
W = addWell(W, G, rock, floor(nx/4)   + floor(1/4*ny)*nx, wellopt{:});
W = addWell(W, G, rock, floor(3*nx/4) + floor(1/4*ny)*nx, wellopt{:});
% production well in the center
wellopt = {'type', 'bhp', 'val', pRef, 'Sign', -1, 'comp_i', [0, 1]};
W = addWell(W, G, rock, floor(nx/2)   + floor(1/4*ny)*nx, wellopt{:});

wellcells = zeros(G.cells.num, 1);
for i = 1 : 3
    wellcells(W(i).cells) = 1;
end
figure
plotCellData(G, wellcells);

facilityModel = FacilityModel(model.fluidModel);
facilityModel = facilityModel.setupWells(W);
model.FacilityModel = facilityModel;
model = model.validateModel(); % setup consistent fields for model (in
                               % particular the facility model for the fluid
                               % submodel)

%% Setup a schedule

clear schedule
schedule.step.val     = [1*day*ones(1, 1); 5*day*ones(2, 1)];
schedule.step.control = (1 : numel(schedule.step.val))';
ratemax = 1*meter^2/day;
currenttime = cumsum(schedule.step.val);
totaltime   = currenttime(end);
for i = 1 : numel(schedule.step.control)
    newrate = ratemax*currenttime(i)/totaltime;
    W(1).val = newrate;
    W(2).val = newrate;
    schedule.control(i) = struct('W', W);
end


%% Run the schedule

[wellSols, states] = simulateScheduleAD(initState, model, schedule);

figure(1); clf
plotToolbar(G, states);


%% Objective functions

% We set up the objective function as beeing the final uplift in one point above
% the production well.

topcell = floor(nx/2) + nx*(ny - 1);
topface = G.cells.faces(G.cells.facePos(topcell) : (G.cells.facePos(topcell + ...
                                                  1)  - 1), :);
topface = topface(topface(:, 2) == 4, 1);
topnode = G.faces.nodes(G.faces.nodePos(topface)); % takes one node from the
                                                   % top face (the first listed)

objUpliftFunc = @(tstep) objUplift(model, states, schedule, topnode, 'tStep', ...
                                     tstep, 'computePartials', true);

%% Compute gradients using the adjoint formulation

adjointGradient = computeGradientAdjointAD(initState, states, model, schedule, objUpliftFunc);


%% There is a routine to check the code using finite difference

% should be used with a smaller schedule otherwise the computation will be
% very long.
compute_numerical_derivative = false;
if compute_numerical_derivative
    objUpliftFunc2 = @(wellSols, states, schedule) computeUplift(model, states, ...
                                                      schedule, topnode, ...
                                                      'computePartials', false);

    fdgrad = computeGradientPerturbationAD(initState, model, schedule, objUpliftFunc2, ...
                                           perturbation);
end

% laststep = numel(states);
% uplift = computeUplift(model, states, schedule, topnode, 'tStep', laststep);
% uplift = uplift{1};