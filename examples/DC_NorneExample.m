%% SPE Example: Simulation of dual-perm poromechancis on a subset of the Norne grid.
% We make use of the fixed-stress split for solution of the coupled problem

%% Load required modules
mrstModule add dual-continuum-mech vemmech ad-core ad-mechanics ad-props

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
perm_fracture = 1000*milli*darcy*ones(G.cartDims); 
rock_fracture = struct('perm', reshape(perm_fracture, [], 1), ...
              'poro', ones(G.cells.num, 1)*0.01);

% matrix
perm_matrix = 0.01*milli*darcy*ones(G.cartDims); 
rock_matrix = struct('perm', reshape(perm_matrix, [], 1), ...
              'poro', ones(G.cells.num, 1)*0.1);
         
          
%% Setup fluid model
fluid = initSimpleADIFluid('phases', 'W', 'mu', 5*centi*poise, 'rho', ...
                                   850*kilogram/meter^3, 'c', ...
                                   1e-8, 'pRef', 0);
fluid_fracture = fluid;
fluid_matrix = fluid; 


%% Setup material parameters for Biot and mechanics
E = 1e9; % Young's module
nu = 0.2; % Poisson's ratio
K_m =  2E9; % Matrix stiffness
K_s = inf; % Solid stiffness

E = repmat(E, G.cells.num, 1);
nu = repmat(nu, G.cells.num, 1);
K_m = repmat(K_m, G.cells.num, 1);
K_s = repmat(K_s, G.cells.num, 1);


%% Shape factor as per  Aziz and Lim, 1995, for use in transfer function
d1 = 0.5; % spacing of fracture set 1 (approx 100 fracs per cell in each direction)
d2 = d1; % spacing of fracture set 2
d3 = d1; % spacing of fracture set 3
a = (pi^2)*(1/(d1^2) + 1/(d2^2) + 1/(d3^2)); % shape factor


%% Setup boundary conditions for mechanics
% We first want to create a structure 'bc', which we can fudge by
% initialising the bc's using pside. 
oside = {'WEST', 'EAST', 'SOUTH', 'NORTH', 'ZMIN', 'ZMAX'};
bc = cell(6,1);
for i = 1:numel(oside)
    bc{i} = pside([], G, oside{i}, 0);
    bc{i} = rmfield(bc{i}, 'type'); 
    bc{i} = rmfield(bc{i}, 'sat');    
end

% DISPLACEMENT BCs
% Find the nodes for the different sides and set the boundaray conditions for
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

% Consolidation problem (zero lateral displacements)
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


% Collect the displacement boundary conditions
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

% FORCE BCs
facesf = bc{6}.face;
force = 50e6; %
n = bsxfun(@times, G.faces.normals(facesf, :), -1./G.faces.areas(facesf));
force_bc = struct('faces', facesf, 'force', bsxfun(@times, n, force));
el_bc = struct('disp_bc', disp_bc, 'force_bc', force_bc);


%% Gravity
% The gravity in this option affects only the fluid behavior
gravity off;
    

%% Setup load for mechanics
% In this example we do not impose any volumetric force
load = @(x) (0*x);


 %% Gather all the mechanical parameters in a struct
mech = struct('E', E, 'nu', nu, 'K_m', K_m, 'K_s', K_s, 'el_bc', el_bc, 'load', load);


%% Fixed stress splitting (FSS)
% fullycoupledOptions = {'verbose', opt.verbose};
% model = DualPermMechWaterModel(G, rock_fracture, fluid_fracture, rock_matrix, fluid_matrix, mech, fullycoupledOptions{: });
% model.transfer_model_object = SimpleTransferFunction();
% model.transfer_model_object.shape_factor_object.shape_factor_value = a;
splittingOptions = {'splittingTolerance', opt.splittingTolerance, ...
                        'splittingVerbose', opt.splittingVerbose};
model = DualPermMechFluidFixedStressSplitModel(G, rock_fracture, fluid_fracture, rock_matrix,...
                                       fluid_matrix, mech, 'fluidModelType', 'water',...
                                       splittingOptions{: });
model.fluidModel.transfer_model_object = SimpleTransferFunction();
model.fluidModel.transfer_model_object.shape_factor_object.shape_factor_value = a;
                                   
%% Setup initial state and fluid BCs
pressure = ones(G.cells.num,1)*0;
state0 = struct('pressure', pressure, 'pom', pressure, 's', ones(G.cells.num, 1), 'swm', ones(G.cells.num, 1));
state0.wellSol = initWellSolAD([], model, state0);
state0.xd = zeros(nnz(~model.mechModel.operators.isdirdofs), 1);
state0 = addDerivedQuantities(model.mechModel, state0);

% Need to initiate the fluid bc's, bc's are the same for micro and macro scales 
bc_f = fluxside([], G, 'WEST', 0, 'sat', 1);
bc_f = fluxside(bc_f, G, 'EAST', 0,'sat', 1);
bc_f = fluxside(bc_f, G, 'SOUTH', 0, 'sat', 1);
bc_f = fluxside(bc_f, G, 'NORTH', 0, 'sat', 1);
bc_f = fluxside(bc_f, G, 'LOWER', 0, 'sat', 1);
bc_f = fluxside(bc_f, G, 'UPPER', 0, 'sat', 1);


%% Setup states and initialise model 
% use empty facilities to when inducing a pressure
model.FacilityModel = FacilityModel(model.fluidModel);
dt = 0.01;
initState = NonLinearSolver().solveTimestep(state0, dt, model, 'bc', bc_f);

% Empty facilities model since we have no wells, meep
W = [];
refdepth = G.cells.centroids(111,3);
prodcell = 61;
W = addWell(W, G, rock_fracture, prodcell, ...
                'Type'    , 'bhp', ...
                'Val'     , 10E6, ...
                'Sign'    , -1,  ...k
                'Comp_i'   , 1, ... % produce water
                'Name'    , 'prod_frac',  ...
                'refDepth', refdepth);

W = addWell(W, G, rock_matrix, prodcell, ...
                'Type'    , 'bhp', ...
                'Val'     , 10E6, ...
                'Sign'    , -1,  ...k
                'Comp_i'   , 1, ... % produce water
                'Name'    , 'prod_mat',  ...
                'refDepth', refdepth);
             
 facilityModel = FacilityModel(model.fluidModel);
 facilityModel = facilityModel.setupWells(W);            
 model.FacilityModel = facilityModel;


%% Setup Schedule and simulate
dt = [1*day*ones(1,5),5*day*ones(1,19)];
schedule = simpleSchedule(dt, 'W', W);

[wellSols, states, schedulereport] = simulateScheduleAD(initState, model, schedule);
states = [{initState};states];
figure
plotToolbar(model.G, states, 'outline', true);
plotWell(G,W, 'fontsize', 0);
axis off
caxis([1E7 3E7])

view(1,40);