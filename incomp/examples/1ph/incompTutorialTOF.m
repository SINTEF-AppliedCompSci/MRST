%% Time-of-flight
% To study how heterogeneity affects  flow patterns and define natural
% time-lines in the reservoir, it is common to study the so-called
% time-of-flight (TOF), i.e., the time it takes an imaginary particle
% released at an inflow boundary or at a perforation of an injector to
% reach a given point in the reservoir. Time-of-flight are usually
% associated with streamline methods, but can also be computed from linear
% steady-state transport equations on the form,
%
% $$v\cdot \nabla \tau = \phi. $$
%
% In this example, we will show how to compute time-of-flight from a given
% flux field.
mrstModule add incomp diagnostics

%% Setup model
% As our model, we will use a logically Cartesian grid in which the
% node perturbed so that grid lines form curves and not lines. For
% simplicity, we assume unit permeability and porosity and use a set of
% source and sink terms to emulate a quater five-point setup.
G = cartGrid([50,50]);
G = twister(G);
G = computeGeometry(G);
rock = makeRock(G, 1, 1);
fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);
src = addSource([], 1, sum(poreVolume(G,rock)), 'sat', 1);
src = addSource(src, G.cells.num, -sum(poreVolume(G, rock)), 'sat', 1);

%% Solve pressure equation
% Compute transmissibilities and solve pressure equation
trans  = computeTrans(G, rock);
xr = incompTPFA(initResSol(G, 0), G, trans, fluid, 'src', src);
clf,
plotCellData(G,xr.pressure, 'edgecolor','k','edgealpha',.05);
title('pressure')
axis equal tight;colormap jet
colorbar

%% Compute time-of-flight
% Once the fluxes are given, the time-of-flight equation is discretized
% with a standard upwind finite-volume method
t0 = tic;
T  = computeTimeOfFlight(xr, G, rock, 'src',src);
toc(t0)

clf,
plotCellData(G, T, 'edgecolor','k','edgealpha',0.05);
title('time-of-flight');
caxis([0,0.8]);axis equal tight;colormap jet
colorbar

%% Copyright notice

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
