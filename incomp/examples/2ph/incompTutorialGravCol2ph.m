%% Two-Phase Flow in Inclined Gravity Column
% In this example, we simulate the injection of a light fluid (CO2) into
% a heavier fluid (brine) inside an inclined sandbox.

mrstModule add incomp


%% Set up model
% To get an inclined reservoir, we manipulate the gravity direction. Since
% gravity is a persistent and global variable, it is important that we
% reset the gravity at the end of the script
rotdeg = 40;
height = 40;
n      = [ 20,  2,  100];
box_sz = [100, 10, 200];
G      = cartGrid(n, box_sz); CG = cartGrid([1 1 1],box_sz);
G      = computeGeometry(G);
rock   = makeRock(G, 0.1*darcy, 0.3);
T      = computeTrans(G, rock, 'verbose', true);
           
fluid  = initSimpleFluid('mu' , [  0.307,   0.049] .* centi*poise     , ...
                         'rho', [973    , 617    ] .* kilogram/meter^3, ...
                         'n'  , [  2    ,   2    ]);

% Redefine gravity direction
rot   = @(theta) makehgtform('xrotate',  theta(1), ...
                             'yrotate', -theta(2), ...
                             'zrotate', -theta(3));
mul   = @(a,b,n) a(1:n,1:n) * reshape(b(1:n), [], 1);
angle = [0, pi*rotdeg/180, 0];
gravity reset on
gravity(mul(rot(angle), gravity(), 3));

%% Compute initial pressure distribution
% Put region of CO2 at bottom of reservoir.
xr = initResSol(G, 1*barsa, 1);
d  = gravity() ./ norm(gravity);
dc = G.cells.centroids * d.';
xr.s(dc>max(dc)-height) = 0;
xr = incompTPFA(xr, G, T, fluid);

cla reset
h = plotCellData(G, convertTo(xr.pressure(1:G.cells.num), barsa), ...
                 'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'FaceAlpha', 0.625);
rotate(h,[0 1 0],rotdeg);
view([0,0]), axis equal tight off
colorbar, colormap(jet(10))

%% Prepare to plot during simulation
clf
h = plotGrid(CG, 'FaceColor', 'none', 'EdgeColor', 'k','LineWidth',1);
rotate(h,[0 1 0],rotdeg);
view([0,0]), axis tight off
hs = plotCellData(G, xr.s, xr.s < .995, 'EdgeColor', 'none');
rotate(hs,[0 1 0],rotdeg);
s  = linspace(0, 1, 64).';
cm = [1-s.^(13/16), 1-s.^6, s.^6];
caxis([0 1]); colormap(cm), colorbar, axis equal

%% Run simulation
% For accuracy, the time step is gradually ramped up
dT = [.5, .5, 1, 1, 1, 2, 2, 2, 5, 5, 10, 10, 15, 20, repmat(25,[1,97])].*day();
t = 0;
for k = 1 : numel(dT),
   xr = implicitTransport(xr, G, dT(k), rock, fluid, 'Verbose', false);

   % Check for inconsistent saturations
   assert (max(xr.s) < 1+eps && min(xr.s) > -eps);

   % Increase time and plot saturation
   t = t + dT(k);
   delete(hs)
   hs = plotCellData(G, xr.s, xr.s <.995, 'EdgeColor', 'none');
   rotate(hs,[0 1 0],rotdeg);
   title(sprintf('%.2f days', t/day));
   drawnow

   % Compute new flow field.
   xr = incompTPFA(xr, G, T, fluid);
end

%% NB: RESET GRAVITY
% Gravity is defined as a persistent and global variable and we therefore
% need to reset it to avoid messing with other examples
gravity reset on

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
