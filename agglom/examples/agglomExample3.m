%% Example 3: Different Flow Indicators
% We present two examples of coarsening algorithms that fall into the
% general algorithmic framework described by Hauge et al [1]:
%
% * uniform refinement of high-flow zones in a Cartesian grid
% * the nononuniform coarsening algorithm
%
% Both algorithms are applied to a 2D five-spot example with heterogeneity
% sampled from a lateral layer of the SPE10 model. For each algorithm, we
% consider three different flow indicators based on the absolute
% permeability, the velocity magnitude, and time-of-flight.
%
% *References:*
%
% # V. L. Hauge, K.-A. Lie, J. R. Natvig,
%   Grid coarsening based on amalgamation for multi-fidelity transport
%   solvers, September 2010,
%   http://www.sintef.no/Projectweb/GeoScale/Publications/

mrstModule add agglom coarsegrid spe10 mimetic incomp diagnostics;

%% Set up and solve flow problem
% As our example, we consider a standard five spot with heterogeneity
% sampled from Model 2 of the 10th SPE Comparative Solution Project.
[G, W, rock] = getSPE10setup(25);
% Set wells to single component mode
for i = 1:numel(W)
    W(i).compi = 1;
end
rock.poro = max(rock.poro, 1e-4);
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
rS = initState(G, W, 0);
S  = computeMimeticIP(G, rock);
rS = incompMimetic(rS, G, S, fluid, 'wells', W);

%%
% We will use the following cell-wise indicators:
%
% * the absolute permeability
% * the magnitude of the velocity at the cell center
% * the product of the forward and backward time-of-flight
%
% For all indicators, we perform a logarithmic scaling
clf
iK = log10(rock.perm(:,1)); iK = iK - min(iK) + 1;
subplot(1,3,1); plotCellData(G,iK,'EdgeColor','none'); axis tight off
title('Permeability')

v  = faceFlux2cellVelocity(G, rS.flux); v  = sqrt(sum(v .^ 2, 2));
iV = log10(v); iV = iV - min(iV) + 1;
subplot(1,3,2); plotCellData(G,iV,'EdgeColor','none'); axis tight off;
title('Velocity')

T  = computeTimeOfFlight(rS, G, rock, 'wells', W);
Tr = computeTimeOfFlight(rS, G, rock, 'wells', W, 'reverse', true);
iT = -log10(T.*Tr); iT = iT - min(iT) + 1;
subplot(1,3,3); plotCellData(G,iT,'EdgeColor','none'); axis tight off;
title('Time-of-flight')

%% Uniform Refinement of High-Flow Zones
% We start by a uniform 3x11 partition of the model and then add an extra
% 3x3 refinement in all blocks in which the block indicator exceeds the
% upper bound: Ib > NU*mean(Ic), where Ib is the indicator per block, Ic
% is the indicator per cell, and NU is a user-prescribed parameter.
p1 = partitionUI(G, [3, 11, 1]);
pK = refineUniformShape(p1, G, iK, 380, 'CartDims', [3,3,1]);
pV = refineUniformShape(p1, G, iV, 390, 'CartDims', [3,3,1]);
pT = refineUniformShape(p1, G, iT, 380, 'CartDims', [3,3,1]);

subplot(1,3,1)
outlineCoarseGrid(G, pK); title(sprintf('Permeability: %d', max(pK)))
subplot(1,3,2)
outlineCoarseGrid(G, pV); title(sprintf('Velocity: %d', max(pV)))
subplot(1,3,3)
outlineCoarseGrid(G, pT); title(sprintf('Time-of-flight: %d', max(pT)))

%%
% Using the two flow-based indicators |iV| and |iT|, rather than the a priori
% permeability indicator |iK|, gives grids that better adapt to the flow
% pattern. This is no surprise, since |iK| only distinguishes between zones
% with high permeability and low permeability, whereas the other two
% indicators also account for the influence of the forces that drive the
% flow.

%% The 'Nonuniform Coarsening' (NUC) Algorithm
% Next, we consider <simpleNUC.html the nonuniform coarsening algorithm>
% and compare the grids generated using the three flow indicators
% introduced above
clf
numBins  = 10;
NL       = 30;
NU       = 80;

I = [iK iV iT];
for i=1:3
   p = segmentIndicator(G, I(:,i), numBins);
   p = mergeBlocks2(p, G, rock.poro, I(:,i), NL, NU);
   p = refineGreedy2(p, G, I(:,i), NU);
   p = mergeBlocks2(p, G, rock.poro, I(:,i), NL, NU);
   subplot(1,3,i), plotCellData(G, I(:,i), 'EdgeColor', 'none')
   outlineCoarseGrid(G,p);
   axis tight off, title(sprintf('%d blocks', max(p)));
end

%%
% In the algorithm above, the greedy algorithm considers a 9-point
% neighbourhood when growing cells. In the original paper, the authors used
% a 5-point neighbourhood, which gives more diamond-shaped cells as shown
% in the figure below
clf
for i=1:2
   p = segmentIndicator(G, iV, numBins);
   p = mergeBlocks(p, G, rock.poro, iV, NL);
   p = refineGreedy2(p, G, iV, NU, 'nlevel',i);
   p = mergeBlocks(p, G, rock.poro, iV, NL);
   subplot(1,2,i), plotCellData(G, iV, 'EdgeColor', 'none')
   outlineCoarseGrid(G,p);
   axis tight off, title(sprintf('%d-neighbour: %d blocks', 2^(i+1)+1, max(p)));
end
%%
% In the two previous plots, we used the greedy algorithm to refine blocks.
% We can, of course, also other types of refinement algorithms of each
% block that exceeds the upper bound on the total flow. For a Cartesian
% grid, the 'refineUniform' method performs a uniform refinement of the
% bounding box of each block that exceeds the upper limit.
clf
for i=1:3
   p = segmentIndicator(G, I(:,i), numBins);
   p = mergeBlocks2(p, G, rock.poro, I(:,i), NL, NU);
   p = refineUniformShape(p, G, I(:,i), NU, 'CartDims',[2,2,1]);
   p = mergeBlocks2(p, G, rock.poro, I(:,i), NL, NU);
   subplot(1,3,i), plotCellData(G, I(:,i), 'EdgeColor', 'none')
   outlineCoarseGrid(G,p);
   axis tight off, title(sprintf('%d blocks', max(p)));
end

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
