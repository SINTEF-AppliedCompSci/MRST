mrstModule add ad-core ad-props ad-blackoil vemmech ad-mechanics

% OPTIONS : They are set by opt, default values are given below.
%
% method_case:
%               'no split'
%               'drained split'
%               'stress split'
%
% norne_case
%               'full'
%               'mini Norne'

% bc_cases
%               'no displacement'
%
%


opt = struct('method_case', 'no split', ...
             'norne_case', 'mini Norne', ...
             'bc_cases', 'no displacement');


%% Load Norne grid
if ~ (makeNorneSubsetAvailable() && makeNorneGRDECL()),
    error('Unable to obtain simulation model subset');
end

is_grid_constructed = false;
if ~is_grid_constructed
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
    %act_org=grdecl.ACTNUM
    grdecl.ACTNUM = ones(size(grdecl.ACTNUM));
    G = processGRDECL(grdecl);
    G = G(1);
    G = computeGeometry(G);
end


%% Setup material parameters from deck for fluid
perm = [grdecl.PERMX, grdecl.PERMY, grdecl.PERMZ];
rock.perm = perm(G.cells.indexMap, :);
rock.poro = max(grdecl.PORO(G.cells.indexMap), 0.1);


is_fluid_initiated = false;
if ~is_fluid_initiated
    pth = getDatasetPath('spe1');
    fn  = fullfile(pth, 'BENCH_SPE1.DATA');
    mrstModule add deckformat
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
end


%% Setup material parameters for Biot and mechanics

E     = 1 * giga * Pascal; % Young's module
nu    = 0.3;               % Poisson's ratio
alpha = 1;                 % Biot's coefficient
Ev = repmat(E, G.cells.num, 1);
nuv = repmat(nu, G.cells.num, 1);
rock.alpha = repmat(alpha, G.cells.num, 1);


%% Setup boundary conditions for mechanics (no displacement)

switch opt.bc_cases
  
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
loadfun = @(x) (0*x);

%% Gravity is on - Only for fluid.
gravity on; 
g = gravity;
if norm(g) == 0
    warning('Gravity off!');
end

%% Setup mechanic struct
mech = struct('Ev', Ev, 'nuv', nuv, 'el_bc', el_bc, 'load', loadfun);

%% Setup model

is_model_setup = false;
if ~is_model_setup
    switch opt.method_case
      case  'no split'
        model = MechBlackOilModel(G, rock, fluid, mech);
      case 'drained split'
        model = MechFluidDrainedSplitModel(G, rock, fluid, mech, 'fluidModelType', 'blackoil');
      case 'stress split'
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, 'fluidModelType', 'blackoil');
      otherwise
        error('method_case not recognized.')
    end
end


%% Setup wells
W = [];
refdepth = G.cells.centroids(1, 3); % for example...
injcell = 10; % for example
W = addWell(W, G, rock, injcell, ...
            'Type'    , 'rate', ...
            'Val'     , 2.5e6/day, ...
            'Sign'    , 1,  ...
            'Comp_i'   , [0, 0, 1], ... % inject gas
            'Name'    , 'inj',  ...
            'refDepth', refdepth);

prodcell = G.cells.num; % for example
W = addWell(W, G, rock, prodcell, ...
            'Type'    ,'bhp', ...
            'Val'     , pRef, ...
            'Sign'    , -1,  ...
            'Comp_i'   , [0, 1, 0], ... % one-phase test case
            'Name'    , 'prod',  ...
            'refDepth', refdepth);

model = model.validateModel();
model.FacilityModel = model.FacilityModel.setupWells(W);


%% Setup schedule
schedule.step.val = [1*day*ones(1, 1); 10*day*ones(11, 1)];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('W', W);

%% Setup initial state
clear initState;
initState.pressure = pRef*ones(G.cells.num, 1);
initState.s = ones(G.cells.num, 1)*[0, 1, 0];
initState.rs = 0.5*fluid.rsSat(initState.pressure);
initState = initDisplacement(model, initState, []);

%% Start simulation

[wsols, states] = simulateScheduleAD(initState, model, schedule);

state = states{end};

figure
plotCellData(G, state.stress(:, 1));
title('First stress component');

figure
plotCellData(G, state.stress(:, 2));
title('Second stress component');

figure
plotCellData(G, state.stress(:, 3));
title('Third stress component');

