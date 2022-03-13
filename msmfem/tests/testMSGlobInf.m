%% Test of use of global information in multiscale basis functions.
% Compare the fine-grid and the multiscale pressure solver by solving the
% single-phase pressure equation
%
% $$\nabla\cdot v = q, \qquad v=-\frac{K}{\mu}\nabla p,$$
%
% for a Cartesian grid with isotropic, homogeneous permeability

try
   require mimetic coarsegrid
catch
   mrstModule add mimetic coarsegrid
end

nx = 20;
ny = 20;
nz = 1;

cellDims  = [nx, ny, nz];
verbose   = mrstVerbose;
G         = cartGrid(cellDims, cellDims);
G         = computeGeometry(G);
rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]);
fluid     = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                            'rho', [1000, 700]*kilogram/meter^3, ...
                            'n'  , [   2,   2]);

% Set two wells, one vertical and one horizontal
W = struct([]);
W = verticalWell(W, G, rock, nx, ny, nz, ...
                 'InnerProduct', 'ip_simple', ...
                 'Type', 'rate', 'Val', 1*meter^3/day, ...
                 'Radius', .1, 'Name', 'I');
W = addWell(W, G, rock, 1, 'Type', 'bhp', 'InnerProduct', 'ip_simple', ...
            'Val', 0, 'Radius', .1, 'Dir', 'x', 'Name', 'P');

src = [];
%src = addSource(src, [1, G.cells.num], [1 -1]./ day());
bc  = [];
%bc = addBc(bc, [1 11], 'flux', [0.5 -0.5]./day());


%% Set up solution structures
% Here we need four solution structures, two for each simulator to hold the
% solutions on the grid and in the wells, respectively.
xRef = initState(G, W, 0);
xMs  = initState(G, W, 0);
xMs2 = initState(G, W, 0);

%% Partition the grid
% We partition the fine grid into a regular 5-by-5-by-2 coarse grid in
% index space so that each coarse block holds 8-by-8-by-5 fine cells. The
% resulting vector <p> has one entry per fine-grid cell giving the index of
% the corresponding coarse block. After the grid is partitioned in index
% space, we postprocess it to make sure that all blocks consist of a
% connected set of fine cells. This step is superfluous for Cartesian
% grids, but is required for grids that are only logically Cartesian (e.g.,
% corner-point and other mapped grids that may contain inactive or
% degenerate cells).
p  = partitionUI(G, [4, 4, 1]);
p  = processPartition  (G, p);

CG = generateCoarseGrid(G, p, 'Verbose', verbose);

S  = computeMimeticIP(G, rock, 'Verbose', verbose);

xRef = incompMimetic(xRef, G, S, fluid, 'wells', W, 'bc', bc, ...
                     'src', src, 'Solver', 'hybrid');

CS  = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]), ...
                           'global_inf', xRef.flux,                ...
                           'Verbose', verbose,'src', src, 'bc', bc);

CS2 = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]), ...
                           'Verbose', verbose, 'src', src,'bc', bc);

W = generateCoarseWellSystem(G, S, CG, CS, ones([G.cells.num, 1]), rock, W);

xMs = solveIncompFlowMS(xMs, G, CG, p, S, CS, fluid, 'wells', W, ...
                        'bc', bc, 'src', src, 'Solver', 'hybrid');

xMs2 = solveIncompFlowMS(xMs2, G, CG, p, S, CS2, fluid, 'wells', W, ...
                         'bc', bc, 'src', src, 'Solver', 'hybrid');

dp = @(x) num2str(convertTo(x.wellSol(1).pressure, barsa));
disp(['DeltaP - Fine:     ', dp(xRef)]);
disp(['DeltaP - Ms_global:', dp(xMs )]);
disp(['DeltaP - Ms_old:   ', dp(xMs2)]);

e   = @(x1, x2) num2str(norm(x1 - x2) ./ norm(x1));
err = @(x1, x2, fld) e(x1.(fld), x2.(fld));
disp(['rel_norm_flux(ref-Ms_global): ', err(xRef, xMs , 'flux')]);
disp(['rel_norm_flux(ref-Ms_old):    ', err(xRef, xMs2, 'flux')]);
disp(['rel_norm_pres(ref-Ms_global): ', err(xRef, xMs , 'pressure')]);
disp(['rel_norm_pres(ref-Ms_old):    ', err(xRef, xMs2, 'pressure')]);

cellNo    = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
plot_var  = @(x) plotCellData(G, x);
plot_flux = @(x) plot_var(accumarray(cellNo, ...
   abs(convertTo(faceFlux2cellFlux(G, x.flux), meter^3/day))));

%% plot output
f = figure;
subplot(3,2,1)
   plot_flux(xRef);
   view(3), camproj perspective, axis tight equal, camlight headlight
   title('Reference solution')
   cax2 = caxis; colorbar

subplot(3,2,3)
   plot_flux(xMs);
   caxis(cax2); colorbar
   title('MS with global info')
   view(3), camproj perspective, axis tight equal, camlight headlight

subplot(3,2,5)
   plot_flux(xMs2);
   caxis(cax2); colorbar
   title('MS without global info')
   view(3), camproj perspective, axis tight equal, camlight headlight

subplot(3,2,2)
   plotCellData(G, convertTo(xRef.pressure(1:G.cells.num), barsa));
   title('Pressure Fine')
   view(3), camproj perspective, axis tight equal, camlight headlight
   cax = caxis; colorbar

subplot(3,2,4)
   plotCellData(G, convertTo(xMs.pressure(1:G.cells.num), barsa));
   title('Pressure Coarse')
   view(3), camproj perspective, axis tight equal, camlight headlight
   caxis(cax); colorbar

subplot(3,2,6)
   plotCellData(G, convertTo(xMs2.pressure(1:G.cells.num), barsa));
   title('Pressure Coarse')
   view(3), camproj perspective, axis tight equal, camlight headlight
   caxis(cax); colorbar

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
