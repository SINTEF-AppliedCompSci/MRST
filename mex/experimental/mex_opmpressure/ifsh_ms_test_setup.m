G = computeGeometry(processGRDECL(makeModel3([100, 60, 15])));
K = logNormLayers(G.cartDims, [10, 300, 40, 0.1, 100]);
rock.perm = bsxfun(@times, [1, 100, 0.1], K(:));
rock.perm = convertFrom(rock.perm(G.cells.indexMap, :), ...
                        milli*darcy);
W = [];
W = verticalWell(W, G, rock, 12, 15, (1:15), ...
                 'Type', 'bhp', 'Val', 500*barsa, 'Radius', 0.125);

W = verticalWell(W, G, rock, 28, 37, (1:15), ...
                 'Type', 'bhp', 'Val', 500*barsa, 'Radius', 0.125);

W = verticalWell(W, G, rock, 85, 38, (1:15), ...
                 'Type', 'bhp', 'Val', 250*barsa, 'Radius', 0.125);

W = verticalWell(W, G, rock, 70, 15, (1:15), ...
                 'Type', 'bhp', 'Val', 300*barsa, 'Radius', 0.125);

fluid = initSingleFluid('mu', 1, 'rho', 1);

xref = initState(G, W, 0, [0, 1]);

t0 = tic;
S    = computeMimeticIP(G, rock, 'verbose', true);
xref = solveIncompFlow(xref, G, S, fluid, 'wells', W)
toc(t0)

src = addSource([], vertcat(W.cells), vertcat(xref.wellSol.flux));

p = mex_partition_ui(double(G.cells.indexMap), G.cartDims, [10, 10, 5]);
p = mex_partition_process(p, G.faces.neighbors);
p = mex_partition_compress(p);

x = initResSol(G, 0);

save ifsh_ms_setup x G rock p src

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.
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
