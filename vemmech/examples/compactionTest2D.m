%% Compaction test case in 2D
%
%  A constant force is imposed at the top and gravity load. The analytical solution
%  can be computed exactly in this case. We compare it with the numerical
%  solution
%
mrstModule add vemmech
%% Define the fluid and rock parameters and set up the grid.

opt = struct('grid_type', 'cartgrid', ... % see squareGrid
            'disturb'   , 0.0, ...      % parameter for disturbing grid
            'E'         , 0.3*1e9, ...  % young's modulus
            'nu'        , 0.4, ...
            'islinear'  , false);       % poisson's ratio

opt.cartDims       = [10 10];
opt.L              = [15*10/10 15];
opt.hanging        = false;
opt.free_top       = true; % If true, the nodes at the top can move freely (no
                           % boundary condition given there)
opt.triangulate    = true; % If true, the horizontal faces are triangulated.
opt.gravity_load   = true; % Use gravity load
opt.top_load       = false;
opt.gtol           = 0.1e-1; % Grid tolerance parameter (used when calling
                             % processGRDECL, see documentation there)
opt.twist          = true;
opt.disturb        = 0.05;
% Different methods are implemented to compute the loading term, see
% paper [Andersen et al: http://arxiv.org/abs/1606.09508v1].
opt.force_method   = 'cell_force'; % 'dual_grad_type'

G = squareGrid(opt.cartDims, opt.L, 'grid_type', opt.grid_type);

if(opt.triangulate)
    face = unique(G.cells.faces(any(bsxfun(@eq, G.cells.faces(:, 2), [3, 4]), 2), 1));
    G = triangulateFaces(G, face);
    G = sortEdges(G);
end

if (opt.twist)
    G = twister(G, opt.disturb);
end

G = createAugmentedGrid(G);
G = computeGeometry(G);

Ev  = repmat(opt.E, G.cells.num, 1);
nuv = repmat(opt.nu, G.cells.num, 1);
C   = Enu2C(Ev, nuv, G);

figure()
clf,
plotGrid(G, 'Marker', '*')
title('Grid')

%% Setup the loads and the boundary conditions
%

% We use the utility function makeCompactionTest
[el_bc, load] = makeCompactionTest(G, opt);


%% Assemble and solve the system

bbsize = 30000-(G.griddim-2)*20000; % block size for the assembly
uu = VEM_linElast(G, C, el_bc, load, ...
                  'blocksize'   , bbsize, ...
                  'force_method', opt.force_method);

figure(),
clf
plotNodeDataDeformed(G, uu(:, G.griddim), uu);
title('Displacement in the vertical direction')

%% Computation of the analytical solution

if(~isempty(el_bc.force_bc))
    ff = abs(el_bc.force_bc.force(1, G.griddim));
else
    ff = 0;
end
start = max(G.faces.centroids(:, G.griddim));
top   = min(G.faces.centroids(:, G.griddim));

fac = 100*300/2;

ana = @(z) ff*(z-start)./(C(1, 1))-double(opt.gravity_load)*fac*((top-start).^2 - (z-top).^2)/C(1, 1);
divana = @(z) (ff./C(1, 1))-double(opt.gravity_load)*fac*(-2*(z-top))/C(1, 1);

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
