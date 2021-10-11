%% Compaction test case 3D
%
%  A constant force is imposed at the top and gravity load. The analytical solution
%  can be computed exactly in this case. We compare it with the numerical
%  solution
%

%% Load required modules

mrstModule add vemmech

%% Define the fluid and rock parameters and set up the grid.

opt = struct('grid_type' , 'triangle', ...
             'disturb'   , 0.0, ...     % parameter for  grid distortion
             'E'         , 0.3*1e9, ... % Young's modulus
             'nu'        , 0.4, ...     % Poison's ratio
             'islinear'  , false);
opt.L            = [15 15 3];
opt.islinear     = false; % if true, the boundary condition is a given linear
                          % displacement, see function makeCompactionTest.
% Different methods are implemented to compute the loading term, see
% paper [Andersen et al: http://arxiv.org/abs/1606.09508v1].
opt.force_method = 'dual_grad_type';
opt.hanging      = false; % If true, zero displacement on the vertical sides.
opt.free_top     = true;  % If true, the nodes at the top can move freely (no
                          % boundary condition given there)
opt.triangulate  = false;  % If true, the horizontal faces are triangulated.
opt.vertical     = false; % Only relevant for norne test case (straightens up
                          % the pillars, see paper [Andersen et al: http://arxiv.org/abs/1606.09508v1])
opt.gravity_load = true;  % Use gravity load
opt.top_load     = true;  % Use force applied at the top
opt.gtol         = 0.1e-1; % Grid tolerance parameter (used when calling
                           % processGRDECL, see documentation there)
opt.ref          = 10;     % Refinement parameter, only used for Norne
opt.flipgrid     = false;  % Rotate the grid (z->x, x->y, y->z) (see paper [Andersen et al: http://arxiv.org/abs/1606.09508v1])
if isempty(mfilename)
    % Running example interactively
    grid_case_number = input(['Choose a grid (type corresponding number): box [1], ' ...
                              'sbed [2], Norne [3]\n']);
else
    % Running example as a function
    grid_case_number = 1;
end
switch grid_case_number
  case 1
    grid_case = 'box';
    opt.cartDims     = [[1 1]*3 10]; % set the Cartesian dimension for the box case
  case 2
    grid_case = 'sbed';
  case 3
    grid_case = 'norne';
  otherwise
    error('Choose grid case by typing number between 1 and 3.');
end

G = complex3DGrid(opt, grid_case);

if (opt.flipgrid)
    G = flipGrid(G);
end
G = createAugmentedGrid(G);
G = computeGeometry(G);

Ev     = repmat(opt.E, G.cells.num, 1);
nuv    = repmat(opt.nu, G.cells.num, 1);
C      = Enu2C(Ev, nuv, G);

figure()
clf;
plotGrid(G);
view(3);

%% Setup the loads and the boundary conditions
if(strcmp(grid_case,'norne'))
    %only rolling in vertical direction this is need since norne has
    %irregular sides and the code do not have genneral implementation of
    %rolling condition at this point
    [el_bc, load] = makeCompactionTest(G, opt, 'rolling_vertical', true);
else
    [el_bc, load] = makeCompactionTest(G, opt);
end


%% Assemble and solve the system

bbsize = 30000-(G.griddim-2)*20000;
lsolve = @mldivide;
fprintf('Running ... ');
uu = VEM_linElast(G, C, el_bc, load, 'linsolve', lsolve, 'blocksize', bbsize, ...
                  'force_method', opt.force_method);
fprintf('done!\n');

figure();
clf;
plotNodeDataDeformed(G, uu(:, 3), uu);
view(3);
title('Vertical displacement')


%% Compute  the analytical solution

ff = abs(el_bc.force_bc.force(1, 3));
start = max(G.faces.centroids(:, 3));
top = min(G.faces.centroids(:, 3));
[lambda, mu] = ENu2LMu_3D(opt.E, opt.nu);
gfac = 10*3000/2; % gravity is 10, density is 3000, 2 is because of derivative
ana = @(z) ff*(z-start)./(C(1, 1))-double(opt.gravity_load)*gfac*((top-start).^2 - (z-top).^2)/C(1, 1);
divana = @(z) (ff./C(1, 1))-double(opt.gravity_load)*gfac*(-2*(z-top))/C(1, 1);

%% Comparison plots

z = G.nodes.coords(:, G.griddim);
z(abs(ana(z))<max(abs(ana(z)))*1e-2) = nan;
zl = unique(z);

figure(),
plot(z, uu(:, G.griddim), '*', zl, ana(zl))
title('Displacement in the vertical direction')
legend({'computed solution', 'analytical solution'})

figure()
zc = G.cells.centroids(:, G.griddim);
div = VEM_div(G);
plot(zc, div*reshape(uu', [], 1)./G.cells.volumes, '*', zc, divana(zc));
title('Divergence');
legend({'computed solution', 'analytical solution'})

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
