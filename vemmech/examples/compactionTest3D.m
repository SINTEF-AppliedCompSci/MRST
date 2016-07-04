%% Compaction test case 3D
% 
% Boundary conditions: rolling condition on left-hand side
%
%
% Possibility to include a top load.
%
% Homogeneous medium, analytical solution, irregular grid

%% Load required modules

mrstModule add vemmech


%% Define parameters

opt = struct('grid_type' , 'triangle', ...
             'disturb'   , 0.0, ...     % parameter for  grid distortion
             'E'         , 0.3*1e9, ... % Young's modulus
             'nu'        , 0.4, ...     % Poison's ratio
             'islinear'  , false);
opt.L            = [15 15 3];
opt.islinear     = false;
opt.force_method = 'dual_grad_type';
opt.hanging      = false;
opt.free_top     = true;  % If true, the nodes at the top can move freely (no
                          % boundary condition given there)
opt.triangulate  = true;  % Triangulate some faces
opt.vertical     = false; % Only relevant for norne test case (straightens up
                          % the pillars, see paper )
opt.gravity_load = true;  % Use gravity load
opt.top_load     = true;  % Use force applied at the top
opt.gtol         = 0.1e-1; % Grid tolerance parameter (used when calling
                           % processGRDECL, see documentation there)
opt.ref          = 10;     % Refinement parameter, only used for Norne
opt.flipgrid     = false;  % Rotate the grid (z->x, x->y, y->z) (see paper )

grid_case_number = input(['Choose a grid (type corresponding number): box [1], ' ...
                    'sbed [2], Norne [3]\n']);
switch grid_case_number
  case 1
    grid_case = 'box';
    opt.cartDims     = [[1 1]*3 10]; % set the Cartesian dimension for the box case
  case 2
    grid_case = 'sbed';
  case 3
    grid_case = 'Norne';
  otherwise
    error('Choose grid case by typing number between 1 and 3.');
end


%% Construct grid

G = complex3DGrid(opt, grid_case);
if (opt.flipgrid)
    G = flipGrid(G);
end
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);

figure()
clf;
plotGrid(G);


%% Setup loads

[el_bc, load] = makeCompactionTest(G, opt);


%% Define rock parameters
Ev     = repmat(opt.E, G.cells.num, 1);
nuv    = repmat(opt.nu, G.cells.num, 1);
C      = Enu2C(Ev, nuv, G);

%% Assemble and solve the system

bbsize = 30000-(G.griddim-2)*20000;
lsolve = @mldivide;
uu = VEM_linElast(G, C, el_bc, load, 'linsolve', lsolve, 'blocksize', bbsize, ...
                  'force_method', opt.force_method);

%% Assemble divergence operator, see paper

div = VEM_div(G);

%% Plotting

figure();
clf;
plotNodeDataDeformed(G, uu(:, 3), uu);cb = colorbar;view(3)
cdiv = div*reshape(uu', [], 1);
axis off

%% Compute  the analytical solution

ff = abs(el_bc.force_bc.force(1, 3));
start = max(G.faces.centroids(:, 3));
top = min(G.faces.centroids(:, 3));
[lambda, mu] = ENu2LMu_3D(opt.E, opt.nu);
ana = @(z) ff*(z-start)./(C(1, 1))+double(opt.gravity_load)*10*300*((z).^2-(start).^2)/C(1, 1);
ana = @(z) ff*(z-start)./(C(1, 1))-double(opt.gravity_load)*50*300*((top-start).^2 - (z-top).^2)/C(1, 1);
divana = @(z) (ff./C(1, 1))-double(opt.gravity_load)*50*300*(-2*(z-top))/C(1, 1);

%% Comparison plots

z = G.nodes.coords(:, 3);
z(abs(ana(z))<max(abs(ana(z)))*1e-2) = nan;
zl = unique(z);
figure(),
subplot(4, 1, 1)
plot(z, uu(:, 3), '*', zl, ana(zl))
subplot(4, 1, 2)
err = (uu(:, 3)-ana(z))./abs(ana(z));
plot(z, err, '*')
subplot(4, 1, 3)
div = VEM_div(G);
plot(G.cells.centroids(:, 3), div*reshape(uu', [], 1)./G.cells.volumes, '*');
subplot(4, 1, 4)
div = VEM_div(G);
zc = G.cells.centroids(:, 3);
diverr = (div*reshape(uu', [], 1)./G.cells.volumes-divana(zc))./abs(divana(zc));
plot(zc, diverr, '*');
