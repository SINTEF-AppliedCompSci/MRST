%% Illustration of the MsRSB method
% This example illustrates the computation of basis functions for the
% multiscale restriction-smoothed basis (MsRSB) method. 
mrstModule add msrsb incomp coarsegrid spe10

%% Set up model
% The mesh is designed to have k*k coarse blocks, that each consists of
% (2m+1)*(2m+1) cells with petrophysical properties sampled from the 15th
% layer of Model 2 of the 10th SPE Comparative Solution Project. By setting
% the flag 'MatrixOutput' to true, we get the linear system as output in
% the fields state.A and state.rhs. We consider a 5x5 coarse mesh with 9x9
% fine cells per coarse block. For simplicity, we only show the basis
% functions in the middle 3x3 coarse blocks to avoid the effect of the
% boundaries.
[k,m] = deal(5,4);
n = k*(2*m+1);
G      = computeGeometry(cartGrid([n n]));
rock   = getSPE10rock(1:n, 1:n, 15);
rock.perm = rock.perm(:, 1);
hT     = computeTrans(G, rock);
bc     = pside([], G, 'West', 100*barsa);
bc     = pside(bc, G, 'East',  50*barsa);
fluid  = initSingleFluid('rho', 1, 'mu', 1);
state0 = initResSol(G, 0);
state  = incompTPFA(state0, G, hT, fluid, 'MatrixOutput', true, 'bc', bc);

%% Create coarse grid and augment with support region
% For the partition, we use a standard load-balanced partition in index
% space. Then we construct a coarse grid structure with appropriate
% geometry information and augment it with support regions
p  = partitionUI(G, [k,k]);
CG = coarsenGeometry(generateCoarseGrid(G,p));
CG = storeInteractionRegion(CG);

%% Predefine matrices used in the iteration
% We construct the matrix that computes one Jacobi increment as well as the
% matrix we postmultiply to normalize the basis functions

% Jacobi increment matrix
[n,m] = deal(n*n, k*k);
J     = spdiags(1./diag(state.A), 0, n, n)*state.A;
w     = 2/3;

% Normalization
lens   = cellfun(@numel, CG.cells.interaction);
blocks = rldecode((1:CG.cells.num)', lens);
ia     = vertcat(CG.cells.interaction{:});
M      = sparse(ia, blocks, ones(size(ia)), CG.parent.cells.num, CG.cells.num);

%% Compute and visualize basis functions
% We show all code lines necessary to compute the basis functions, perform
% 100 iterations, and visualize how the basis function for the coarse block
% at the midpoint of the domain evolves over time.

% Set initial state for basis functions and plot the basis in coarse block
% number 13, i.e., the middle block of the domain.
P = zeros(n, m);
for i=1:m, P(p==i,i)=1; end
subplot(2,2,1);
plotCellData(G,P(:,13),P(:,13)>0,'EdgeColor','none'); 
plotFaces(CG,1:CG.faces.num,'EdgeColor','w','LineWidth',1); 
axis equal tight off, axis([9 36 9 36])

% Run 100 steps of the iteration and show basis function in block 13 after
% step 3, 10, and 100.
plt = [0 3 10 100 nan];
pp = 2;
for j=1:100
    incr = J*P;
    incr = incr.*M;
    P    = P - w*incr;
    P    = bsxfun(@rdivide, P, sum(P, 2));
    if j==plt(pp)
        subplot(2,2,pp);
        plotCellData(G,P(:,13),P(:,13)>0,'EdgeColor','none'); 
        plotFaces(CG,1:CG.faces.num,'EdgeColor','w','LineWidth',1); 
        axis equal tight off
        axis([9 36 9 36])
        pp=pp+1;
    end
end

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
