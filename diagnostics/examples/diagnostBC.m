%% Illustrate time-of-flight/tracer for a case with flow across boundaries
mrstModule add diagnostics incomp 


%% Set up and solve flow problem
% We use a square reservoir with one well and pressure BC on one side

% Grid
[nx,ny] = deal(100);
G = cartGrid([nx,ny,1],[250,250,10]);
G = computeGeometry(G);

% Petrophysical data
p = gaussianField(G.cartDims(1:2), [0.01 0.4], [11 3], 150);
K = p.^3.*(1.5e-5)^2./(0.81*72*(1-p).^2);

rock = makeRock(G, K(:), p(:));

hT  = computeTrans(G, rock);
pv  = sum(poreVolume(G,rock));

% Fluid model
gravity reset off
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

% Wells
W = addWell([],  G, rock, round(nx/3)+nx*round(ny/3), 'Type', 'rate', 'Comp_i', 1, 'name', 'I1', 'Val', 10/day);

% Pressure BC
bc  = pside([], G, 'Right', 100.*barsa());

% Initial reservoir state
state = initState(G, W, 0.0, 1.0);

%% Compute basic quantities
state = incompTPFA(state, G, hT, fluid, 'wells', W, 'bc', bc);

% split flow boundary into sub-regions (10 faces in each)
bpart = ceil((1:numel(bc.face))'/10);
% compute TOF/tracer (setting computeWellTOFs to true will also compute 
% boundary region TOFs)

D = computeTOFandTracer(state, G, rock, 'wells', W, 'bc', bc, ...
                        'partitionBoundary', bpart, ...
                        'computeWellTOFs', true);

figure, 
subplot(1,2,1)
plotCellData(G, D.outpart, 'EdgeColor', 'none');
title('Boundary tracer region partition'), axis tight
subplot(1,2,2)
plotCellData(G, D.outtracer(:,9), 'EdgeColor', 'none')
title('Seventh backward tracer region '), axis tight

%% Copyright notice

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
