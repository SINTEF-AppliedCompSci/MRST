% simpleMSBCExample
% The example shows how to put Dirichlet boundary conditions on the system.
mrstModule add coarsegrid mimetic msmfem incomp


cellDims = [40, 40, 10];
verbose  = true;

gravity off;

G = cartGrid(cellDims, cellDims);
G = computeGeometry(G);

rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(0.3            , [G.cells.num, 1]);
W = struct([]);
W = verticalWell(W, G, rock, 40, 40, 1:10, ...
                 'InnerProduct', 'ip_simple', ...
                 'Type', 'rate', 'Val', 1*meter^3/day, ...
                 'Radius', .1, 'Name', 'I');
W = addWell(W, G, rock, 1:40, 'Type','bhp', ...
            'InnerProduct', 'ip_simple', ...
            'Val', 0, 'Radius', .1, 'Dir', 'x', 'Name', 'P');

fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1000, 700]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);
xrRef = initResSol(G, 0.0);

p  = partitionUI(G, [5, 5, 2]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p, 'Verbose', verbose);
S  = computeMimeticIP(G, rock, 'Verbose', verbose);
%
bc = pside([], G, 'LEFT', 0);
CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]), ...
                          'Verbose', verbose, 'bc', bc);

W = generateCoarseWellSystem(G, S, CG, CS, ones([G.cells.num, 1]), rock, W);

xRef = initState(G, W, 0);
xMs  = initState(G, W, 0);

xRef = incompMimetic    (xRef, G, S, fluid, 'bc', bc, 'wells', W);
xMs  = solveIncompFlowMS(xMs, G, CG, p, S, CS, fluid, 'wells', W, ...
                         'bc', bc);

%% plot output
f = figure;
cellNo = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';
plot_var = @(x) plotCellData(G, x, 'EdgeAlpha', 0.125);
plot_pres = @(x) plot_var(convertTo(x.pressure(1:G.cells.num), barsa));
plot_flux = @(x) plot_var(log10(accumarray(cellNo, ...
   abs(convertTo(faceFlux2cellFlux(G, x.flux), meter^3/day)))));
subplot(2,2,1)
   plot_pres(xRef); title('Pressure Fine [bar]')
   view(3), camproj perspective, axis tight equal, camlight headlight
   cax = caxis; colorbar

subplot(2,2,2)
   plot_pres(xMs); title('Pressure Coarse [bar]')
   view(3), camproj perspective, axis tight equal, camlight headlight
   caxis(cax);
   colorbar

subplot(2,2,3)
   plot_flux(xRef); title('Flux intensity Fine')
   view(3), camproj perspective, axis tight equal, camlight headlight
   cax2 = caxis; colorbar

subplot(2,2,4)
   plot_flux(xMs); title('Flux intensity Coarse')
   view(3), camproj perspective, axis tight equal, camlight headlight
   caxis(cax2); colorbar

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
