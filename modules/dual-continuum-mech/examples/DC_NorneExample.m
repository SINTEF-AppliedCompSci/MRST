%%  Example of a poromechanical dual-continuum simulation on a 3D geological 
%   grid with void fractures
%   
%   The geological grid is from the Norne field (https://opm-project.org/?page_id=559)
%   with a top boundary stress and free vertical displacements at all boundaries,
%   apart from the bottom, which stays fixed. Drainage comes from a centrally located
%   well, into which we allow flow from both the matrix and fracture continua. 
%   Finally, we consider the void-space constitutive coefficient model (see 
%   example_void_fractures).
%
%   As a note, in this example, following equilibration, it is fairly
%   difficult to see a dual-continuum response with the given boundary conditions
%   (well setup). In this case, one could justify a single-porosity model, 
%   however, for different BCs this may not be the case. This example
%   therefore serves as a test bed for investigating when we should and
%   shouldn't consider DC responses at 3D geological scales. 
%

%% Load required modules
mrstModule add dual-continuum-mech dual-porosity vemmech ad-core ad-mechanics ad-props mrst-gui

%% Setup default options
clear opt
opt = struct( 'fluid_model'        , 'water', ...
              'method'             , 'fixed stress splitting', ...
              'nonlinearTolerance' , 1e-6, ...
              'splittingTolerance' , 1e-6, ...
              'verbose'            , true, ...
              'splittingVerbose'   , true);


%% Setup Grid - We use a subset of the Norne grid, which is freely available
% via the open porous media project: https://opm-project.org/ 
grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
fullCartDims = grdecl.cartDims;
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);
grdecl = cutGrdecl(grdecl, [10 20; 35 45; 1 5]);

grdecl.ACTNUM = ones(size(grdecl.ACTNUM));
G = processGRDECL(grdecl);
G = G(1);
G = computeGeometry(G);   

plotGrid(G)
view(3)


%% Setup rock parameters (for flow)
%fracture
perm_fracture = 10*darcy*ones(G.cartDims); 
rock_fracture = struct('perm', reshape(perm_fracture, [], 1), ...
              'poro', ones(G.cells.num, 1)*0.01,...
              'vol_fraction', ones(G.cells.num, 1)*0.01);

% matrix
perm_matrix = 0.01*milli*darcy*ones(G.cartDims); 
rock_matrix = struct('perm', reshape(perm_matrix, [], 1), ...
              'poro', ones(G.cells.num, 1)*0.1,...
              'vol_fraction', ones(G.cells.num, 1)-rock_fracture.vol_fraction);
         
          
%% Setup fluid model
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                                   1000*kilogram/meter^3, 'c', ...
                                   4e-10, 'pRef', 0);
fluid_fracture = fluid;
fluid_matrix = fluid; 


%% Setup material parameters for Biot and mechanics
E = 1E8; % Young's module
nu = 0.2; % Poisson's ratio
E_m =  1E9; % Matrix stiffness
nu_m = 0.2;
K_s = inf; % Solid stiffness

E = repmat(E, G.cells.num, 1);
nu = repmat(nu, G.cells.num, 1);
E_m = repmat(E_m, G.cells.num, 1);
nu_m = repmat(nu_m, G.cells.num, 1);
K_s = repmat(K_s, G.cells.num, 1);


%% Setup boundary conditions for mechanics
% we first want to create a structure 'bc', which we can fudge by
% initialising the bc's using pside. 
oside = {'WEST', 'EAST', 'SOUTH', 'NORTH', 'ZMIN', 'ZMAX'};
bc = cell(6,1);
for i = 1:numel(oside)
    bc{i} = pside([], G, oside{i}, 0);
    bc{i} = rmfield(bc{i}, 'type'); 
    bc{i} = rmfield(bc{i}, 'sat');    
end

% DISPLACEMENT BCs
% Find the nodes for the different sides and set the boundary conditions for
% elasticity.
for i = 1 : 6
    inodes = mcolon(G.faces.nodePos(bc{i}.face), G.faces.nodePos(bc{i}.face + 1) - 1);
    nodes = unique(G.faces.nodes(inodes));
    disp_bc = struct('nodes'   , nodes,      ...
                     'uu'      , 0,          ...
                     'faces'   , bc{i}.face, ...
                     'uu_face' , 0,          ...          
                     'mask'    , true(numel(nodes), G.griddim));
    bc{i}.el_bc = struct('disp_bc', disp_bc, 'force_bc', []);
end
bcdisp_zero = @(x) x*0.0; % Boundary displacement function set to zero.

% consolidation problem (zero lateral displacements)
bc_el_sides{1} = bc{1}; 
bc_el_sides{1}.el_bc.disp_bc.mask(:, 3) = false;   % x fixed, y fixed, z free
bc_el_sides{2} = bc{2}; 
bc_el_sides{2}.el_bc.disp_bc.mask(:, 3) = false;   % x fixed, y fixed, z free
bc_el_sides{3} = bc{3}; 
bc_el_sides{3}.el_bc.disp_bc.mask(:, 3) = false;   % x fixed, y fixed, z free
bc_el_sides{4} = bc{4}; 
bc_el_sides{4}.el_bc.disp_bc.mask(:, 3) = false;   % x fixed, y fixed, z free
bc_el_sides{5} = bc{5}; 
bc_el_sides{5}.el_bc.disp_bc.mask(:, :) = true;   % x fixed, y fixed, z fixed
bc_el_sides{6} = bc{6}; 
bc_el_sides{6}.el_bc.disp_bc.mask(:, 3) = false;   % x fixed, y fixed, z free

% collect the displacement boundary conditions
nodes = [];
faces = [];
mask = [];
for i = 1 : numel(bc)
    if(~isempty(bc_el_sides{i}))
        nodes = [nodes; bc_el_sides{i}.el_bc.disp_bc.nodes]; %#ok
        faces = [faces; bc_el_sides{i}.el_bc.disp_bc.faces]; %#ok
        mask  = [mask; bc_el_sides{i}.el_bc.disp_bc.mask]; %#ok
    end
end

disp_node = bcdisp_zero(G.nodes.coords(nodes, :));
disp_faces = bcdisp_zero(G.faces.centroids(faces, :)); 
disp_bc = struct('nodes', nodes, 'uu', disp_node, 'faces', faces, 'uu_face', disp_faces, 'mask', mask); 

% force BCs
facesf = bc{6}.face;
force = 20E6; %
n = bsxfun(@times, G.faces.normals(facesf, :), -1./G.faces.areas(facesf));
force_bc = struct('faces', facesf, 'force', bsxfun(@times, n, force));
el_bc = struct('disp_bc', disp_bc, 'force_bc', force_bc);


%% Gravity
% the gravity in this option affects only the fluid behaviour
gravity off;
    

%% Setup load for mechanics
% in this example we do not impose any volumetric force
load = @(x) (0*x);


%% Gather all the mechanical parameters in a struct
mech = struct('E', E, 'nu', nu, 'E_m', E_m, 'nu_m', nu_m, 'K_s', K_s, 'el_bc', el_bc, 'load', load);


%% Setup model object
fullycoupledOptions = {'verbose', true};
DC_model = DualContMechWaterModel(G, {rock_fracture, rock_matrix}, {fluid_fracture, fluid_matrix}, mech, fullycoupledOptions{:});
d1 = 1; % spacing of fracture set 1
d2 = d1; % spacing of fracture set 2
d3 = d2;
fracture_spacing = repmat([d1,d2,d3],G.cells.num,1);
shape_factor_name = 'Lim_AzizShapeFactor';
DC_model.transfer_model_object = SimpleTransferFunction(shape_factor_name, fracture_spacing);
                                   

%% Initialise model
pressure = zeros(G.cells.num,1);
state0 = struct('pressure', pressure, 'pressure_matrix', pressure,...
                's', ones(G.cells.num, 1), 'swm', ones(G.cells.num, 1));
state0.xd = zeros(nnz(~DC_model.mechModel.operators.isdirdofs), 1);
state0.wellSol = initWellSolAD([], DC_model, state0);
state0 = addDerivedQuantities(DC_model.mechModel, state0);

% need to initiate the fluid bc's, bc's are the same for micro and macro scales 
bc_f = fluxside([], G, 'WEST', 0, 'sat', 1);
bc_f = fluxside(bc_f, G, 'EAST', 0,'sat', 1);
bc_f = fluxside(bc_f, G, 'SOUTH', 0, 'sat', 1);
bc_f = fluxside(bc_f, G, 'NORTH', 0, 'sat', 1);
bc_f = fluxside(bc_f, G, 'LOWER', 0, 'sat', 1);
bc_f = fluxside(bc_f, G, 'UPPER', 0, 'sat', 1);

% initialise in response to applied load at the top boundary
DC_model = DC_model.validateModel(); % instantiate fluidModel.FacilityModel
DC_model.FacilityModel = FacilityModel(DC_model.fluidModel); % map facility model from fluidModel up to main model
initState = NonLinearSolver().solveTimestep(state0, 0.01, DC_model, 'bc', bc_f);

% allow matrix and fracture continua to equilibrate
N = 1;
while norm(initState.pressure - initState.pressure_matrix, 1)>0.1
    initState = NonLinearSolver().solveTimestep(initState, 10*day, DC_model, 'bc', bc_f);
    N = N + 1;
    if N > 100
        break
    end
end


%% Setup Schedule and simulate
% wells
refdepth = G.cells.centroids(111,3);
prodcell = 61;
W = addWell([], G, rock_fracture, prodcell, ...
                'Type'    , 'bhp', ...
                'Val'     ,  0, ...
                'Sign'    , -1,  ...k
                'Comp_i'   , 1, ... % produce water
                'Name'    , 'prod_frac',  ...
                'refDepth', refdepth);
% facility model object is processed by parent class but assigned to current
% model object
DC_model.FacilityModel = FacilityModel(DC_model.fluidModel).setupWells(W);

% setup schedule and simulate
dt = [1*day*ones(1,5),5*day*ones(1,30)];
schedule = simpleSchedule(dt, 'bc', bc_f, 'W', W);
[wellSols, states, schedulereport] = simulateScheduleAD(initState, DC_model, schedule);
states = [{initState};states];

%% visualise
figure
plotToolbar(DC_model.G, states, 'outline', true);
plotWell(G,W, 'fontsize', 0);
axis off
caxis([0 1.8E7])
view(1,40);

plotWellSols(wellSols,'field','qWr')