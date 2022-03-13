%% MsFV example on PEBI grid
% This example demonstrates the use of the MsFV method on unstructured
% grids using a PEBI grid. This grid is unstructured in ij-space.
%
% We start by loading the required modules and a stored binary file that
% includes the grids. If the stored binary file is not found, it is
% downloaded from the Internet (approximately 15 MB download)

mrstModule add msfvm mrst-gui mex agglom coarsegrid incomp

if ~exist('pebi_with_dual.mat', 'file')
    disp('Did not find grid, downloading (~15.0 MB)')
    msfvmdir = fullfile(mrstPath('query', 'msfvm'), 'examples');
    unzip('http://www.sintef.no/project/MRST/pebi_with_dual.mat.zip', msfvmdir)
end

load pebi_with_dual
%% Plot the coarse grids (dual and primal)
v = [150, 70];

clf; plotCellData(G, mod(CG.partition,17), 'edgecolor', 'k', 'edgea', .5);
title('Primal partition')
view(v), axis tight off;

figure;
plotCellData(G, mod(CG.partition,17), DG.lineedge, 'EdgeColor', 'w')
title('Dual grid edges')
view(v), axis tight off;
%% Set up boundary conditions and fluid
bc = [];
d  = abs(G.faces.centroids(:,1) - max(G.faces.centroids(:,1)));
es = find (d < 1e6*eps);

d  = abs(G.faces.centroids(:,1) - min(G.faces.centroids(:,1)));
ws = find (d < 1e6*eps);

bc = addBC(bc, ws, 'pressure', 1000*barsa());
bc = addBC(bc, es, 'pressure', 250*barsa());
fluid        = initSingleFluid('mu' ,    1*centi*poise     , ...
                               'rho', 1014*kilogram/meter^3);
clf; plotGrid(G, 'facea', 0)
plotFaces(G, ws, 'FaceColor', 'red')
plotFaces(G, es, 'FaceColor', 'blue')
title('Boundary conditions (red=1000 bar, blue=250 bar)')
view(v)
axis tight off

%% Set up the permeability
% The values here can be changed to solve either a homogenous permeability
% field or two variants of the SPE10 mapped onto the grid using nearest
% neighbor interpolation.
%
% The default is the Tarbert layers which is inbetween the homogenous
% permeability and the Ness layers in terms of difficulty.

gridtype = 'spe10';
% Uncomment to take uniform permeability instead
% gridtype = 'uniform';

% Take the Tarbert layers 
offset = 0;
% Uncomment to take Ness layers instead
% offset = 35;

switch lower(gridtype)
    case 'uniform'
        rock.perm = 100*milli*darcy*ones(G.cells.num,1);
        rock.poro = ones(G.cells.num, 1);
    case 'spe10'
        % We have 12 layers in the model
        Nz = 12;
        mrstModule add spe10 mex

        [Gspe, W, rockspe] = getSPE10setup((1:Nz) + offset);
        % Normalize coordinates to make it easier to map the values
        normalize = @(x) bsxfun(@rdivide, x, max(x));
        c = normalize(Gspe.cells.centroids);
        rock.perm = zeros(G.cells.num, 3);
        for i = 1:3
            % Interpolate the permeability using nearest neighbor
            % interpolation
            pperm = TriScatteredInterp(c, rockspe.perm(:, i) , 'nearest');
            rock.perm(:, i) = pperm(normalize(G.cells.centroids));
        end
    otherwise
        error('Unknown case')
end

sol = initState(G, [], 0, 1);
T = computeTrans(G, rock);

% Plot the permeability
figure(1);
clf;
plotCellData(G, log10(rock.perm(:, 1)));
axis tight off; view(v); colorbar;
title('log10 permx')

%% Solve the reference solution
% Standard two point scheme
disp 'Solving TPFA reference'
tic;
solRef = incompTPFA(sol, G, T, fluid, 'bc', bc, 'MatrixOutput', true);
toc;
%% Set up msfv solver and solve
% We call the solver twice to both obtain the uniterated solution as well
% as a solution where a single smoother cycle has been applied (5 steps
% with block gauss-seidel / DMS).

solvems = @(iterator, iterations, smoother, subiterations, omega, verb, restart) ... 
          solveMSFV_TPFA_Incomp(sol, G, CG, T, fluid, 'Dual', DG, 'Reconstruct', true, 'bc', bc, 'Verbose', verb, 'SpeedUp', true,...
          'Iterator', iterator, 'Iterations', iterations, 'Smoother', smoother, 'Subiterations', subiterations, 'Omega', omega, 'Restart', restart);

disp 'Solving MSFVM'
[sol_msfv] = solvems('msfvm', 1, 'dms', 5, 1, true, 0);



disp 'Solving MSFVM'
[sol_msfv_0] = solvems('msfvm', 0, 'dms', 0, 1, true, 0);

%% Plot all solutions in seperate plots and the error
figure(1);
clf;
plotCellData(G, solRef.pressure);
axis tight off; view(v); colorbar;
set(gca, 'CLim', [0, max(solRef.pressure)]);
title('TPFA')

figure(2);
clf;
plotCellData(G, sol_msfv.pressure);
axis tight off; view(v); colorbar;
    title('MsFVM with one smoother cycle')
% Use same axis scaling as the TPFA solution
set(gca, 'CLim', [0, max(solRef.pressure)]); 

figure(3);
clf;
plotCellData(G, sol_msfv_0.pressure);
axis tight off; view(v); colorbar;
    title('MsFVM initial solution')
% Use same axis scaling as the TPFA solution
set(gca, 'CLim', [0, max(solRef.pressure)]); 


figure(4);
clf;
err = abs(sol_msfv.pressure - solRef.pressure)./solRef.pressure;
plotCellData(G, err);
outlineCoarseGrid(G, CG.partition);
axis tight off; view(v); colorbar;
title('Error and primal coarse grid')

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
