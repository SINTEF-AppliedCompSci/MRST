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
