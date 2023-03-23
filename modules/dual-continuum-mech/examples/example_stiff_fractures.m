%% Example of a poromechanical dual-continuum simulation on a 2D grid with
%  stiff fractures
%
% The test problem is a uniaxial single-phase consolidation problem, such that the
% bottom boundary is fixed, and the left, right and top boundaries admit
% only vertical displacement. All of the boundaries are no flux boundaries,
% apart from the top boundary which is drained. 
%
% Within this example we demonstrate the following solvers:
%   * 'fully coupled'          : fully coupled solver


%% Load required modules
mrstModule add dual-continuum-mech ad-core ad-mechanics dual-porosity ad-props vemmech


%% Setup default options
opt = struct('cartDims'            , [2, 2], ...
             'L'                  , [10, 10], ...
             'fluid_model'        , 'water', ...
             'verbose'            , false);

             
%% Setup Grid
G = cartGrid(opt.cartDims, opt.L);
G = computeGeometry(G);
%plotGrid(G);


%% Setup rock parameters (for flow)
% fracture
perm_fracture = 1000*milli*darcy*ones(G.cartDims); 
rock_fracture = struct('perm', reshape(perm_fracture, [], 1), ...
              'poro', ones(G.cells.num, 1)*0.003,... 
              'vol_fraction', ones(G.cells.num, 1)*0.005); 
               % note, poro = intrinsic_poro*vol_fraction, therefore, for
               % this case intrinsic_poro is 0.6

% matrix
perm_matrix = 0.01*milli*darcy*ones(G.cartDims); 
rock_matrix = struct('perm', reshape(perm_matrix, [], 1), ...
              'poro', ones(G.cells.num, 1)*0.1,...
              'vol_fraction', ones(G.cells.num, 1)-rock_fracture.vol_fraction);
          
          
%% Setup fluid model, the current implementation only admits single-phase flow
% but the code is sufficiently generalisable to admit multi-phase flow. 
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                                   1000*kilogram/meter^3, 'c', ...
                                   4e-10, 'pRef', 0);
fluid_fracture = fluid;
fluid_matrix = fluid; 


%% Setup material parameters for Biot and mechanics
nu = 0.2; % Poisson's ratio
E_m = 36e9; % Young's modulus of matrix continuum
nu_m = nu; % Poisson's ratio of matrix continuum
E_f = 3.6E6; % Young's modulus of matrix continuum
nu_f = nu; % Poisson's ratio of fracture continuum
K_s = 70E9; % Solid stiffness
E = HS_bound(E_m, E_f, nu_m, nu_f, rock_matrix.vol_fraction(1),...
             rock_fracture.vol_fraction(1), 'lower'); % Young's modulus of fractured rock mass with HS lower bound

E = repmat(E, G.cells.num, 1);
nu = repmat(nu, G.cells.num, 1);
E_m = repmat(E_m, G.cells.num, 1);
nu_m = repmat(nu_m, G.cells.num, 1);
E_f = repmat(E_f, G.cells.num, 1);
nu_f = repmat(nu_f, G.cells.num, 1);
K_s = repmat(K_s, G.cells.num, 1);


%% Setup boundary conditions for mechanics
% we first want to create a structure 'bc', which we can fudge by
% initialising the bc's using pside. 
oside = {'WEST', 'EAST', 'SOUTH', 'NORTH'};
bc = cell(4,1);
for i = 1:numel(oside)
    bc{i} = pside([], G, oside{i}, 0);
    bc{i} = rmfield(bc{i}, 'type'); 
    bc{i} = rmfield(bc{i}, 'sat');    
end

% Displacement BCs
% Find the nodes for the different sides and set the boundaray conditions for
% elasticity.
for i = 1 : 4
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
bc_el_sides{1} = bc{1}; 
bc_el_sides{1}.el_bc.disp_bc.mask(:, 2) = false;   % x fixed, y free
bc_el_sides{2} = bc{2}; 
bc_el_sides{2}.el_bc.disp_bc.mask(:, 2) = false;   % x fixed, y free
bc_el_sides{3} = bc{3}; 
bc_el_sides{3}.el_bc.disp_bc.mask(:, :) = true;    % x fixed, y fixed
bc_el_sides{4} = bc{4}; 
bc_el_sides{4}.el_bc.disp_bc.mask(:, 2) = false;   % x fixed, y free

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

% Force BCs
facesf = bc{4}.face;
force = 2e6;
force_bc = struct('faces', facesf, 'force', force*ones(G.cartDims(1),1)*[0 -1]);
el_bc = struct('disp_bc', disp_bc, 'force_bc', force_bc);


%% Gravity
% the gravity in this option affects only the fluid behaviour
gravity off;
    

%% Setup load for mechanics
% in this example we do not impose any volumetric force
load = @(x) (0*x);


%% Gather all the mechanical parameters in a struct that will be used to
% to define the mech problem
mech = struct('E', E, 'nu', nu, 'E_m', E_m, 'nu_m', nu_m, 'E_f', E_f, 'nu_f', nu_f,...
              'K_s', K_s, 'el_bc', el_bc, 'load', load);


%% Setup fully coupled and fixed stress splitting models
fullycoupledOptions = {'verbose', opt.verbose};
DC_model = DualContMechWaterModel(G, {rock_fracture, rock_matrix}, {fluid_fracture, fluid_matrix}, mech, fullycoupledOptions{:});
d1 = 0.1; % spacing of fracture set 1
d2 = d1; % spacing of fracture set 2
fracture_spacing = repmat([d1,d2],G.cells.num,1);
shape_factor_name = 'Lim_AzizShapeFactor';
DC_model.transfer_model_object = SimpleTransferFunction(shape_factor_name, fracture_spacing);
DC_model = DC_model.validateModel();


%% Setup initial state and fluid BCs
pressure = zeros(G.cells.num,1);
state0 = struct('pressure', pressure, 'pressure_matrix', pressure, 's', ones(G.cells.num, 1), 'swm', ones(G.cells.num, 1));

% need to initiate the fluid bc's, bc's are the same for micro and macro scales 
bc_f0 = fluxside([], G, 'WEST', 0, 'sat', 1);
bc_f0 = fluxside(bc_f0, G, 'EAST', 0,'sat', 1);
bc_f0 = fluxside(bc_f0, G, 'SOUTH', 0, 'sat', 1);
bc_f0 = fluxside(bc_f0, G, 'NORTH', 0, 'sat', 1);


%% Simulate 
time = [0, logspace(-4,2,60)];
dt = diff(time);
[p_m, p_f, ~, states] = simDC_mech(state0, dt, DC_model, bc_f0);


%% Plot results
figure
semilogx(time, sum(p_m,1)./G.cells.num, '-', 'linewidth', 1.5)
hold on
semilogx(time, sum(p_f,1)./G.cells.num, '-d', 'linewidth', 1.5, 'markersize', 6)
hold on
xlabel('time [s]')
ylabel('average pressure [Pa]')
legend('matrix', 'fracture')
title('Results for the intrinsic fracture stiffness simulation')

        