%% Example: Coupling between fluid flow (blackoil model) and rock mechanics (elasticity)
% We use the same blackoil equations as in the ad-blackoil model. For rock
% mechanics, we use the VEM discretization that is given in the vemmech module
%

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui


%% Options for this example

opt.norne_case = 'mini Norne'; % 'full' or 'mini Norne'
%
% 'full'       : 7392 cells
% 'mini Norne' :  605 cells

opt.bc_case    = 'bottom fixed'; % 'no displacement' or 'bottom fixed'
%
% 'no displacement' : All nodes belonging to external faces have displacement
%                     equal to zero
% 'bottom fixed'    : The nodes that belong to the bottom have zero
%                     displacement, while a given pressure is imposed on
%                     the external faces that are not bottom faces.

opt.method     = 'fully coupled'; % 'fully coupled' 'fixed stress splitting'
%
% 'fully coupled'          : The mechanical and flow equations are solved fully couplde.
% 
% 'fixed stress splitting' : The mechanical and flow equations are solved
%                            sequentially using a fixed stress splitting

opt.fluid_model = 'single phase'; % 'blackoil' 'single phase'
% 
% 'blackoil'     : blackoil model is used for the fluid (gas is injected, see
%                  schedule below)
% 'single phase' : single phase model is used for the fluid 


%% Load Norne grid

if ~ (makeNorneSubsetAvailable() && makeNorneGRDECL()),
    error('Unable to obtain simulation model subset');
end


grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
fullCartDims = grdecl.cartDims;
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);
switch opt.norne_case
  case 'full'
    grdecl = cutGrdecl(grdecl, [10 25; 35 55; 1 22]);
  case 'mini Norne'
    grdecl = cutGrdecl(grdecl, [10 20; 35 45; 1 5]);
  otherwise
    error('norne_case not recognized');
end
grdecl.ACTNUM = ones(size(grdecl.ACTNUM));
G = processGRDECL(grdecl);
G = G(1);
G = computeGeometry(G);


%% Setup fluid parameters from SPE1

perm = [grdecl.PERMX, grdecl.PERMY, grdecl.PERMZ];
rock.perm = perm(G.cells.indexMap, :);
rock.poro = max(grdecl.PORO(G.cells.indexMap), 0.1);

pth = getDatasetPath('spe1');
fn  = fullfile(pth, 'BENCH_SPE1.DATA');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);
fluid = initDeckADIFluid(deck);
fluid = rmfield(fluid, 'relPerm');
fluid = rmfield(fluid, 'pcOW');
fluid = rmfield(fluid, 'pcOG');

% Setup quadratic relative permeabilities, since SPE1 relperm are a bit rough.
fluid.krW = @(s) s.^2;
fluid.krG = @(s) s.^2;
fluid.krOW = @(s) s.^2;
fluid.krOG = @(s) s.^2;
pRef = deck.PROPS.PVTW(1);


%% Setup material parameters for Biot and mechanics

E          = 1 * giga * Pascal; % Young's module
nu         = 0.3;               % Poisson's ratio
alpha      = 1;                 % Biot's coefficient
Ev         = repmat(E, G.cells.num, 1);
nuv        = repmat(nu, G.cells.num, 1);
rock.alpha = repmat(alpha, G.cells.num, 1);


%% Setup boundary conditions for mechanics (no displacement)

switch opt.bc_case

  case 'no displacement'

    ind = (G.faces.neighbors(:, 1) == 0 | G.faces.neighbors(:, 2) == 0);
    ind = find(ind);
    nodesind = mcolon(G.faces.nodePos(ind), G.faces.nodePos(ind + 1) - 1);
    nodes = G.faces.nodes(nodesind);
    bcnodes = zeros(G.nodes.num);
    bcnodes(nodes) = 1;
    bcnodes = find(bcnodes == 1);
    nn = numel(bcnodes);
    u = zeros(nn, 3);
    m = ones(nn, 3);
    disp_bc = struct('nodes', bcnodes, 'uu', u, 'mask', m);
    force_bc = [];

  case 'bottom fixed'

    nx = G.cartDims(1);
    ny = G.cartDims(2);
    nz = G.cartDims(3);

    % Find the bottom nodes. On these nodes, we impose zero displacement

    c = zeros(nx*ny*nz, 1);
    c(G.cells.indexMap) = (1 : numel(G.cells.indexMap))';
    bottomcells = c(nx*ny*(nz - 1) +  (1 : (nx*ny))');
    faces = G.cells.faces(mcolon(G.cells.facePos(bottomcells), G.cells.facePos(bottomcells ...
                                                      + 1) - 1), :);
    bottomfaces = faces( faces(:, 2) == 6  , 1);
    indbottom_nodes = mcolon(G.faces.nodePos(bottomfaces), G.faces.nodePos(bottomfaces ...
                                                      + 1) - 1);
    bottom_nodes = G.faces.nodes(indbottom_nodes);
    isbottom_node = false(G.nodes.num, 1);
    isbottom_node(bottom_nodes) = true;
    bcnodes = find(isbottom_node);

    nn = numel(bcnodes);
    u = zeros(nn, 3);
    m = ones(nn, 3);
    disp_bc = struct('nodes', bcnodes, 'uu', u, 'mask', m);

    % Find outer faces that are not at the bottom. On these faces, we impose
    % a given pressure.

    is_outerface1 = (G.faces.neighbors(:, 1) == 0);
    is_outerface1(bottomfaces) = false;
    is_outerface2 = G.faces.neighbors(:, 2) == 0;
    is_outerface2(bottomfaces) = false;

    is_outerface = is_outerface1 | is_outerface2;

    outer_faces = find(is_outerface);

    outer_pressure = pRef;
    signcoef = (G.faces.neighbors(outer_faces, 1) == 0) - ...
        (G.faces.neighbors(outer_faces, 2) == 0);
    n = bsxfun(@times, G.faces.normals(outer_faces, :), signcoef./ ...
               G.faces.areas(outer_faces));
    force = bsxfun(@times, n, outer_pressure);

    force_bc = struct('faces', outer_faces, 'force', force);


  otherwise
    error('bc_cases not recognized')
end

el_bc = struct('disp_bc' , disp_bc, ...
               'force_bc', force_bc);


%% Setup load for mechanics

% In this example we do not impose any volumetric force
loadfun = @(x) (0*x);



%% Gather all the mechanical parameters in a struct

mech = struct('Ev', Ev, 'nuv', nuv, 'el_bc', el_bc, 'load', loadfun);


%% Gravity
% The gravity in this option affects only the fluid behavior
gravity on;


%% Setup model

modeltype = [opt.method, ' and ', opt.fluid_model];
switch modeltype
  case 'fully coupled and blackoil'
    model = MechBlackOilModel(G, rock, fluid, mech);
  case 'fixed stress splitting and blackoil'
    model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, 'fluidModelType', ...
                                           'blackoil');
  case 'fully coupled and single phase'
    model = MechSinglephaseModel(G, rock, fluid, mech);
  otherwise
    error('modeltype not recognized.');
end



%% Setup wells
W = [];
refdepth = G.cells.centroids(1, 3); % for example...
injcell  = 10; % for example...
prodcell = G.cells.num; % for example...

W = addWell(W, G, rock, injcell, ...
            'Type'    , 'rate', ...
            'Val'     , 2.5e6/day, ...
            'Sign'    , 1,  ...
            'Comp_i'   , [0, 0, 1], ... % inject gas
            'Name'    , 'inj',  ...
            'refDepth', refdepth);

W = addWell(W, G, rock, prodcell, ...
            'Type'    ,'bhp', ...
            'Val'     , pRef, ...
            'Sign'    , -1,  ...
            'Comp_i'   , [0, 1, 0], ... % one-phase test case
            'Name'    , 'prod',  ...
            'refDepth', refdepth);

switch opt.fluid_model
  case 'blackoil'
    W(1).compi = [0, 0, 1];
    W(2).compi = [0, 1, 0];
  case 'single phase'
    W(1).compi = [1];
    W(2).compi = [1];
  otherwise
    error('fluid_model not recognized.')
end

facilityModel = FacilityModel(model.fluidModel); 
facilityModel = facilityModel.setupWells(W);
model.FacilityModel = facilityModel;


%% Setup schedule
schedule.step.val     = [1*day*ones(1, 1); 10*day*ones(11, 1)];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control      = struct('W', W);

%% Setup initial state
clear initState;
initState.pressure = pRef*ones(G.cells.num, 1);
initState.s        = ones(G.cells.num, 1)*[0, 1, 0];
initState.rs       = 0.5*fluid.rsSat(initState.pressure);
initState          = computeInitDisp(model, initState, []);


%% Start simulation

[wsols, states] = simulateScheduleAD(initState, model, schedule);

figure
plotToolbar(G, states, 'outline', true);
view([7, 40]);
