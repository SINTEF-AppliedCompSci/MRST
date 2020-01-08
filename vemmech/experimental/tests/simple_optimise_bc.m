mrstModule add vemmech


%% defining model parameters
E = 0.3*1e9;
nu = 0.45;
res = 40; % grid resolution
L = 15;
unit_force = 30000;

%% Creating grid
G = cartGrid([res, res], [L, L]);
G = computeGeometry(G);
G = createAugmentedGrid(G);

%% setting material parameters
Ev = repmat(E, G.cells.num, 1);
nuv = repmat(nu, G.cells.num,1);
C = Enu2C(Ev, nuv, G);

cfun = @(u) C;


%% Setup load and boundary condition

% compute top faces
tmp = pside([], G, 'ymax', 1);
top_faces = tmp.face;
tmp = pside([], G, 'xmax', 1);
side_faces = tmp.face;


% compute bottom nodes
tmp = pside([], G, 'ymin', 1);
bottom_faces = tmp.face;
bottom_nodes = unique(G.faces.nodes(mcolon(G.faces.nodePos(bottom_faces), ...
                                           G.faces.nodePos(bottom_faces+1)-1)));
% set zero displacement at bottom
% disp_bc = SparseTensor(zeros(numel(bottom_nodes), 1), bottom_nodes, 'n') * ...
%           SparseTensor([1, 1], 'd');

disp_eps = 1e-4;
disp_gen = @(u) u * (linspace(-disp_eps, disp_eps, res+1))';

% set constant force at top
top_force_bc = SparseTensor(unit_force) * ...
               SparseTensor([], top_faces, 'f') * ...
               SparseTensor([0, 1], 'd');

side_force_bc = SparseTensor(unit_force) * ...
                SparseTensor([], side_faces, 'f') * ...
                SparseTensor([0, 1], 'd');

bcfun = @(u) struct('disp_bc', SparseTensor(disp_gen(u(3)), bottom_nodes, 'n') * ...
                               SparseTensor([1 0], 'd'), ...
                    'force_bc', (SparseTensor(u(1)) * top_force_bc) + ...
                                (SparseTensor(u(2)) * side_force_bc));
% bcfun = @(u) struct('disp_bc', disp_bc, ...
%                     'force_bc', (SparseTensor(u(1)) * top_force_bc) + ...
%                                 (SparseTensor(u(2)) * side_force_bc));

loadfun = @(u) @(x) 0 * x; % no load

%% compute "correct" solution
u_correct = [0.192823, 0.0114290, 0.83]';
%u_correct = [0, 0, 0.83]';
[dd, extra] = VEM_linElast_AD(G, cfun(u_correct), bcfun(u_correct), loadfun(u_correct));
dd = dd';
%dd = dd(~extra.disc.isdirdofs);

%% define objective function
scaling = 1e9;
obj_fun = @(u, x) ...
          match_displacement_objfun(u, x, extra.disc.isdirdofs, dd(:), scaling);

%% Solve optimization problem
%u = 0.10;
u = [0.20, 0.20, 0.20]';

[foptval, uopt, history, uu_opt, extra] = ...
    optimize_mech(u, G, bcfun, cfun, loadfun, obj_fun);

uopt
u_correct
uopt - u_correct

% %% Assemble and solve system
% [uu, extra] = VEM_linElast_AD(G, C, el_bc, load);

% %% Show result
fac = 1e4;
figure;
plotNodeDataDeformed(G, uu_opt(:, G.griddim), fac * uu_opt, 'edgealpha', 0.1);
% title('Displacement in the vertical direction');