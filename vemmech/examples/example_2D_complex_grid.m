%% A short example solving linear elasticity on a complex grids
%

%% Load required modules

mrstModule add vemmech


%% Define parameters
%
opt = struct('L'         , [1 1], ...
             'cartDims'  , [4 4], ...
             'grid_type' , 'square', ...
             'disturb'   , 0.02, ... %parameter for disturbing grid
             'E'         , 4e7, ...  %youngs modolo
             'nu'        , 0.44);% poiso ratio

%% make a mixed type of grid
%

G = squareGrid(opt.cartDims, opt.L, 'grid_type', 'mixed4', 'disturb', opt.disturb);
G = removeCells(G, [140 : 151, 233 : 235, 249 : 250, 93, 117 : 118]');
G = createAugmentedGrid(G);
G = computeGeometry(G);
figure()
plotGrid(G);
title('Original grid');

%% Find sides of domain
%

[Lx, Ly] = deal(opt.L(1), opt.L(2));
assert(G.griddim == 2);
x = [0, Lx];
bc = cell(4,1);
for i = 1 : 2
    faces = find(abs(G.faces.centroids(:, 1) - x(i)) < eps);
    bc{i} = addBC([], faces, 'pressure', 0);
    bc{i} = rmfield(bc{i}, 'type');
    bc{i} = rmfield(bc{i}, 'sat');
end
y = [0, Ly];
for i = 1 : 2
    faces = find(abs(G.faces.centroids(:, 2) - y(i)) < eps);
    bc{i + 2} = addBC([], faces, 'pressure', 0);
    bc{i + 2} = rmfield(bc{i + 2}, 'type');
    bc{i + 2} = rmfield(bc{i + 2}, 'sat');
end

for i = 1 : 4
    inodes = mcolon(G.faces.nodePos(bc{i}.face), G.faces.nodePos(bc{i}.face + 1) - 1);
    nodes = unique(G.faces.nodes(inodes));
    % Set up the boudary conditions.
    % Note that the fields 'faces' and 'uu_face' are not needed for VEM but
    % will become necessary when running other methods (such as MPSA)
    disp_bc = struct('nodes', nodes, ...
        'uu', 0, ...
        'faces', bc{i}.face, ...
        'uu_face', 0, ...
        'mask', true(numel(nodes), G.griddim));
    bc{i}.el_bc = struct('disp_bc', disp_bc, 'force_bc', []);
end

%% Set up the loading term
%
% Gravity is our loading term
density = 3000; %kg/m^3
grav = 10; %gravity
load = @(x) - (grav*density)*repmat([0, 1], size(x, 1), 1);

%% Set up the displacement boundary condtions
%

% Set boundary displacement function to zero
bcdisp = @(x) x*0.0;

% Set up the boundary conditions for each side
bc_el_sides{1} = bc{1}; % x side is fixed
bc_el_sides{2} = bc{2}; % y side is fixed
bc_el_sides{3} = [];    % bottom  is free
bc_el_sides{4} = [];    % top is free

% Collect the displacement boundary conditions
nodes = [];
faces = [];
mask  = [];
for i = 1 : numel(bc)
    if (~isempty(bc_el_sides{i}))
        nodes = [nodes; bc_el_sides{i}.el_bc.disp_bc.nodes]; %#ok
        faces = [faces; bc_el_sides{i}.el_bc.disp_bc.faces]; %#ok
        mask  = [mask; bc_el_sides{i}.el_bc.disp_bc.mask];   %#ok
    end
end
disp_node  = bcdisp(G.nodes.coords(nodes, :));
disp_faces = bcdisp(G.faces.centroids(faces, :));
disp_bc = struct('nodes', nodes, 'uu', disp_node, 'faces', faces, 'uu_face', disp_faces, 'mask', mask);

%% Set up the force boundary conditions
%

% A force is applied on the top surface. It is discontinuous in the sense
% that it takes two different values on the left and on the right.
force = 50*barsa;
face_force = @(x) force*sign(x(:, 1) - opt.L(1)/2) + 100*barsa;
faces = bc{4}.face;
% Set up the force boundary structure,
% Note that the unit for the  force is  Pa/m^3.
force_bc = struct('faces', faces, 'force', bsxfun(@times, face_force(G.faces.centroids(faces, :)), [0 -1]));


% Final structure defining the boundary conditions
el_bc = struct('disp_bc', disp_bc, 'force_bc', force_bc);

%% Define the rock parameters
%
Ev = repmat(opt.E, G.cells.num, 1);
nuv = repmat(opt.nu*0 + 0.4, G.cells.num, 1);
C = Enu2C(Ev, nuv, G);

%% Assemble and solve the system
%
lsolve = @mldivide;
[uu, extra] = VEM_linElast(G, C, el_bc, load, 'linsolve', lsolve);

%% Plot displacement in y direction
%
plotopts = {'EdgeAlpha', 0.0, 'EdgeColor', 'none'};
plotopts1 = {'EdgeAlpha', 0.01};
figure()
plotNodeData(G, uu(:, 2), plotopts{ : });
colorbar();
title('Displacement in the y direction');

%% Plot the displacement in x direction
%
figure()
plotNodeData(G, uu(:, 1), plotopts{ : });
title('Displacement in the x direction');
colorbar();

%% plot the deformed grid
%
figure()
fac = 1;
plotGridDeformed(G, uu*fac); axis tight
title('The deformed grid');

%% Calculate and plot the divergence
%
vdiv = VEM_div(G);
mdiv = vdiv*reshape(uu', [], 1)./G.cells.volumes;
figure()
plotCellDataDeformed(G, mdiv, uu*fac, plotopts1{ : }); axis tight
colorbar()
title('Divergence of the displacement field');

%% Calulate the stress and strain
%
op = extra;
stress = reshape(op.D*op.WC'*op.assemb'*reshape(uu', [], 1), 3, [])';

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
