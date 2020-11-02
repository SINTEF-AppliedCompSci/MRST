%% Example of a poromechanical dual-continuum simulation on a 2D grid with
%  an anisotropic arrangement of the fracture continuum (which is itself
%  isotropic). 
%
% The test problem is a uniaxial single-phase consolidation problem, such that the
% bottom boundary is fixed, and the left, right and top boundaries admit
% only vertical displacement. All of the boundaries are no flux boundaries,
% apart from the top boundary which is drained. 
%
% Within this example we demonstrate the following solvers:
%   * 'fully coupled'          : fully coupled solver 
%


%% Load required modules
mrstModule add dual-continuum-mech vemmech ad-core ad-mechanics ad-props dual-porosity


%% Setup default options
    opt = struct('cartDims'            , [2, 2], ...
                 'L'                  , [10, 10], ...
                 'fluid_model'        , 'water', ...
                 'verbose'            , false);


%% Setup Grid
G = cartGrid(opt.cartDims, opt.L);
G = createAugmentedGrid(G);
G = computeGeometry(G);
plotGrid(G);


%% Setup rock parameters (for flow)
% fracture
perm_fracture = [0, 1000]*milli*darcy.*ones(G.cells.num,G.griddim);
rock_fracture = struct('perm', perm_fracture, ...
              'poro', ones(G.cells.num, 1)*0.002,...
              'vol_fraction', ones(G.cells.num, 1)*0.0025);

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
E_m = 36e9; E_f = 3.6e6; % Continuum Young's moduli
[E1, E2] = Reuss_Voigt_bounds(E_m, E_f, rock_matrix.vol_fraction(1), rock_fracture.vol_fraction(1));
nu_m = 0.2; nu_f = 0.2; % Continuum Poisson's ratio
[nu12, nu21] = Reuss_Voigt_bounds(nu_m, nu_f, rock_matrix.vol_fraction(1), rock_fracture.vol_fraction(1));
mu_m = E_m/(2*(1+nu_m)); mu_f = E_f/(2*(1+nu_f)); % Continuum shear moduli
[mu12, ~] = Reuss_Voigt_bounds(mu_m, mu_f, rock_matrix.vol_fraction(1), rock_fracture.vol_fraction(1));
K_s = 70E9; % Solid stiffness

E = repmat([E1, E2], G.cells.num, 1);
nu = repmat(nu21, G.cells.num, 1);
mu = repmat(mu12, G.cells.num, 1);
E_m = repmat(E_m, G.cells.num, 1);
nu_m = repmat(nu_m, G.cells.num, 1);
E_f = repmat([E_f, E_f], G.cells.num, 1);
nu_f = repmat(nu_f, G.cells.num, 1); 
mu_f = repmat(mu_f, G.cells.num, 1); % ficticious property here
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
% Find the nodes for the different sides and set the boundary conditions for
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
bc_el_sides{3}.el_bc.disp_bc.mask(:, :) = true;   % x fixed, y fixed, z free
bc_el_sides{4} = bc{4}; 
bc_el_sides{4}.el_bc.disp_bc.mask(:, 2) = false;   % x fixed, y freee


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
% The gravity in this option affects only the fluid behaviour
gravity off;
    

%% Setup load for mechanics
% in this example we do not impose any volumetric force
load = @(x) (0*x);


 %% Gather all the mechanical parameters in a struct that will be used to
 % to define the mech problem
mech = struct('E', E, 'nu', nu, 'mu', mu, 'E_m', E_m, 'nu_m', nu_m, 'E_f', E_f, 'nu_f', nu_f,...
              'mu_f', mu_f, 'K_s', K_s, 'el_bc', el_bc, 'load', load);
[mech.C, mech.invC, mech.invCi]= AnisoEnu2C(E, nu, mu, G);


%% Setup model
fullycoupledOptions = {'verbose', opt.verbose};
DC_model = DualContMechWaterModel(G, {rock_fracture, rock_matrix}, {fluid_fracture, fluid_matrix}, mech, fullycoupledOptions{:});
d1 = 0.1; % spacing of fracture set 1
a = (pi^2)*(1/(d1^2)); % shape factor
DC_model.transfer_model_object = SimpleTransferFunction();
DC_model.transfer_model_object.shape_factor_object.shape_factor_value = a;
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
time = logspace(-4,2,60);
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
title('Results for the anisotropic fractured rock mass simulation')

          