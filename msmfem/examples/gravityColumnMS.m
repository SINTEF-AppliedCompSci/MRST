%% Multiscale Pressure Solver: Simple Case Driven by Gravity
% Compare the fine-scale and the multiscale solver for the single-phase
% pressure equation
%
% $$\nabla\cdot v = q, \qquad
%    v=\textbf{--}\frac{K}{\mu} \bigl[\nabla p+\rho g\nabla z\bigr],$$
%
% This example is a direction continuation of <gravityColumn.html "My First
% Flow-Solver: Gravity Column"> and introduces the multiscale flow solver
% without going into specific details. More details can be found in the
% <simpleBCMS.html "Basic Multiscale Tutorial">.

mrstModule add coarsegrid mimetic msmfem incomp
%% Define the model
% The domain [0,1]x[0,1]x[0,30] is discretized using a Cartesian grid with
% homogeneous isotropic permeability of 100 mD. The fluid has density 1000
% kg/m^3 and viscosity 1 cP and the pressure is 100 bar at the top of the
% structure
nx = 2; ny = 2; nz = 30;
Nx = 1; Ny = 1; Nz =  6;
gravity reset on
G          = cartGrid([nx, ny, nz], [1, 1, 30]);
G          = computeGeometry(G);
rock.perm  = repmat(0.1*darcy(), [G.cells.num, 1]);
fluid      = initSingleFluid('mu' ,    1*centi*poise     , ...
                             'rho', 1014*kilogram/meter^3);

bc         = pside([], G, 'TOP', 100.*barsa);

%% Assemble and solve the fine-scale linear system
S   = computeMimeticIP(G, rock);
sol = incompMimetic(initResSol(G , 0.0), G, S, fluid, 'bc', bc);

%% Plot the fine-scale solution
newplot;
subplot(3, 2, [1 3])
   plotCellData(G, convertTo(sol.pressure(1:G.cells.num), barsa), ...
                'EdgeColor', 'k');
   set(gca, 'ZDir', 'reverse'), title('Fine-scale pressure [bar]')
   view(45,5), cx = caxis; colorbar

%% Multiscale system
p  = partitionUI(G, [Nx, Ny, Nz]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p);

CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]),'bc', bc);
xrMs = solveIncompFlowMS (initResSol(G, 0.0), G, CG, p, ...
                          S, CS, fluid, 'bc', bc);

%% Plot the coarse-scale solution
% As one clearly can see from the plot, the multiscale solution only
% captures the gravity effect on the coarse scale. To also capture
% fine-scale gravity effects, one can either add extra correction functions
% or insert the multiscale solution into the fine-scale equations and solve
% for a residual correction
subplot(3, 2, [2 4])
   plotCellData(G, convertTo(xrMs.pressure(1:G.cells.num), barsa), ...
                'EdgeColor', 'k');
   set(gca, 'ZDir', 'reverse'); title('Coarse-scale pressure [bar]')
   view(45,5); caxis(cx); colorbar

subplot(3, 2, [5 6]);
   plot(1:nz, convertTo(sol .pressure(1:nx*ny:nx*ny*nz), barsa()), '-o',...
        1:nz, convertTo(xrMs.pressure(1:nx*ny:nx*ny*nz), barsa()), '-*');
   legend('fine','coarse','Location','NorthWest');

%%
displayEndOfDemoMessage(mfilename)

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
