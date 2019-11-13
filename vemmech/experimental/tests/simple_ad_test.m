mrstModule add vemmech

%% defining model parameters
E = 0.3*1e9;
nu = 0.45;
res = 100; % grid resolution
L = 15;
top_force = -30000;

% setting ADI
top_force = initVariablesADI(top_force);

%% Creating grid
G = cartGrid([res, res], [L, L]);
G = computeGeometry(G);
G = createAugmentedGrid(G);

%% setting material parameters
Ev = repmat(E, G.cells.num, 1);
nuv = repmat(nu, G.cells.num,1);
C = Enu2C(Ev, nuv, G);

%% Setup load and boundary condition

% compute top faces
tmp = pside([], G, 'ymax', 1);
top_faces = tmp.face;

% compute bottom nodes
tmp = pside([], G, 'ymin', 1);
bottom_faces = tmp.face;
bottom_nodes = unique(G.faces.nodes(mcolon(G.faces.nodePos(bottom_faces), ...
                                           G.faces.nodePos(bottom_faces+1)-1)));
% set zero displacement at bottom
el_bc.disp_bc = struct('nodes', bottom_nodes, ...
                       'mask', ones(numel(bottom_nodes), 2), ...
                       'uu', zeros(numel(bottom_nodes), 2));


% set constant force at top
el_bc.force_bc = SparseTensor(top_force) * ...
                 SparseTensor([], top_faces, 'f') * ...
                 SparseTensor([0, 1], 'd');

load = @(x) 0*x; % no load

%% Assemble and solve system
[uu, extra] = VEM_linElast_AD(G, C, el_bc, load);

%% Show result
fac = 1000;
plotNodeDataDeformed(G, uu(:, G.griddim), fac * uu, 'edgealpha', 0.2);
title('Displacement in the vertical direction');