%% Example 4: Nonuniform Coarsening of SPE10
% This example shows how to perform flow-based coarsening using the
% algorithm proposed by Aarnes et al. [1], in which the key idea is to
% distinguish cells in high-flow and low-flow regions. To this end, we
% start by the following ad-hoc algorithm:
%
% # Segment log(|v)| to get an initial partition.
% # Merge blocks whose volumes are below a lower limit.
% # Refine blocks in which the flow exceeds an upper limit
% # Merge blocks whose volumes are below a lower limit.
%
% The original algorithm consisted of the four steps above. Here, however,
% we repeat the last two steps a few times to obtain better grids. The
% resulting algorithm is a special case of a more general framework
% described by Hauge et al. [2].
%
% The purpose of the example is to demonstrate the main coarsening steps
% for a single layer of the SPE10 model and compare the flow-adapted grid
% with streamlines traced from the flow field used to adapt the grid.
%
% *References:*
%
% # J. E. Aarnes, V. L. Hauge, Y. Efendiev, Coarsening of three-dimensional
%   Structured and unstructured grids for subsurface flow. Advances in
%   Water Resources, Volume 30, Issue 11, November 2007, pp. 2177--2193.
% # V. L. Hauge, K.-A. Lie, J. R. Natvig,
%   Grid coarsening based on amalgamation for multi-fidelity transport
%   solvers, September 2010.
%   http://www.sintef.no/Projectweb/GeoScale/Publications/

mrstModule add agglom coarsegrid spe10 incomp;

%% Set up and solve flow problem
% As our example, we consider a standard five spot with heterogeneity
% sampled from Model 2 of the 10th SPE Comparative Solution Project, which
% we assume that the user has downloaded to a specific data directory using
% the functions supplied as part of the MRST data sets.
[G, W, rock] = getSPE10setup(25);
rock.poro    = max(rock.poro, 1e-3);
rock.perm    = rock.perm(:,1);
fluid        = initSingleFluid('mu' , 1*centi*poise     , ...
                               'rho', 1014*kilogram/meter^3);

% Solve for pressure and velocity, using no-flow boundary conditions.
rSol = initState(G, W, 0);
T    = computeTrans(G, rock);
rSol = incompTPFA(rSol, G, T, fluid, 'wells', W);

%% Compute indicator and set parameters used later in the algorithm
% Compute a scalar, cell-wise velocity measure and use the logarithm of the
% velocity measure as our flow indicator
v = faceFlux2cellVelocity(G, rSol.flux);
v = sqrt(sum(v .^ 2, 2));
iVel = log10(v);
iVel = iVel - min(iVel) + 1;
NL = 10;
NU = 75;

%% Segment indicator into bins
% In the first step, we segment the indicator value into ten bins. An
% alternative choice could have been to set |numBins =
% round(max(iVel)-min(iVel))|
p = segmentIndicator(G, iVel, 10);
plotCoarseningStep(p, G, rock.poro, iVel, NL, NU, 1, true);

%% Merge small blocks
% The segmentation will typically create a speckle of small blocks that we
% do not want in our coarse grid. Each small block is therefore merged with
% its neighbouring block that has the closest indicator value. Here, we use
% porosity as indicator when measuring which blocks to merge
p2 = mergeBlocks2(p, G, rock.poro, iVel, NL, NU);
plotCoarseningStep(p2, G, rock.poro, iVel, NL, NU, 2, true);

%% Refine blocks
% Merging blocks may give new blocks through which there is too much flow.
% We therefore refine blocks in which the flow exceeds an upper bound.
p = refineGreedy2(p2, G, iVel, NU);
plotCoarseningStep(p, G, rock.poro, iVel, NL, NU, 3, true);

%% Merge small blocks
% The refinement step may have created blocks that have too small volume.
% As a safeguard, we merge these
p = mergeBlocks2(p, G, rock.poro, iVel, NL, NU);
[b, p] = findConfinedBlocks(G, p);
plotCoarseningStep(p, G, rock.poro, iVel, NL, NU, 4, true);

%% Repeat refinement and merging
p = refineGreedy2(p, G, iVel, NU);
p =  mergeBlocks2(p, G, rock.poro, iVel, NL, NU);
plotCoarseningStep(p, G, rock.poro, iVel, NL, NU, 4, true);

%%
p = refineGreedy2(p, G, iVel, NU);
p = mergeBlocks2(p, G, rock.poro, iVel, NL, NU);
plotCoarseningStep(p, G, rock.poro, iVel, NL, NU, 4, true);

%% Other refinement algorithms
% The |refineGreedy| routine grows blocks somewhat agressively by adding
% rings of neighbouring cells at the time. The results are typically worse
% than for |refineGreedy2|
p = refineGreedy(p2, G, iVel, NU);
[b, p] = findConfinedBlocks(G, p);
p = mergeBlocks2(p, G, rock.poro, iVel, NL, NU);
plotCoarseningStep(p, G, rock.poro, iVel, NL, NU, 4, 1);

%%
% Better results may be obtained if we use the |refineGreedy3| method
% in which the neighbouring cells are sorted in descending order in terms
% of the number of faces shared with cells in the growing block.
% Unfortunately, the method is quite expensive and its use is not
% recommended for large models.
p = refineGreedy3(p2, G, iVel, NU, 'nlevel',2);
[b, p] = findConfinedBlocks(G, p); %#ok<*ASGLU>
p = mergeBlocks2(p, G, rock.poro, iVel, NL, NU);
plotCoarseningStep(p, G, rock.poro, iVel, NL, NU, 4, 1);

%%
% Alternatively, we can perform a Cartesian refinement of each block, which
% leads to a higher number of blocks.
pC = refineUniformShape(p2, G, iVel, NU,'cartDims', [2 2 1]);
pC = mergeBlocks2(pC, G, rock.poro, iVel, NL, NU);
[b, p] = findConfinedBlocks(G, p);
plotCoarseningStep(pC, G, rock.poro, iVel, NL, NU, 4, 1);

%%
% Or, we can perform a refinement in which we try to impose a fixed
% default partition in all parts that exceed the upper bound. And if this
% is not sufficient, we further subdivide the remaining violating blocks
% rectangularly.
pf = partitionUI(G, [6 22 1]);
pC = refineUniformShape(p2, G, iVel, NU,'fixPart', pf, 'CartDims',[2 1 1]);
pC = mergeBlocks2(pC, G, rock.poro, iVel, NL, NU);
[b, p] = findConfinedBlocks(G, p);
plotCoarseningStep(pC, G, rock.poro, iVel, NL, NU, 4, 1);


%% Compare with streamlines
% To show how the NUC grid adapts to the flow field, we plot permeability
% field, pressure field with streamlines, and flow-adapted grid. As we can
% observe from the figure, the flow-adapted grid aligns neatly with the
% high-flow and low-flow regions
clf,
subplot(1,3,1)
plotCellData(G,log10(rock.perm),'EdgeColor','none'); axis tight off

subplot(1,3,2)
mrstModule add streamlines
plotCellData(G,rSol.pressure,'FaceAlpha',.5,'EdgeColor','none');
nx = G.cartDims(1);
cells = [...
   (-nx/6:nx/6) + W(5).cells - 4*nx, (-nx/6:nx/6) + W(5).cells + 4*nx, ...
   (0:nx/10) + W(1).cells + 4*nx, (-nx/10:0)+W(2).cells + 4*nx, ...
   (0:nx/10) + W(3).cells - 4*nx, (-nx/10:0)+W(4).cells - 4*nx, ...
   ];
h = streamline(pollock(G, rSol, cells')); set(h,'Color','b');
rSol.flux = -rSol.flux;
h = streamline(pollock(G, rSol, cells')); set(h,'Color','b');
rSol.flux = -rSol.flux;
plotWell(G,W); axis tight off

subplot(1,3,3)
plotCellData(G,iVel,'EdgeColor','none')
outlineCoarseGrid(G, p); axis tight off

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
