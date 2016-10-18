load 'PebiGrid3D'

tol = 1e-3;
f   = boundaryFaces(G);
top = abs(G.faces.centroids(f,3)-1) < tol;
bot = abs(G.faces.centroids(f,3)  ) < tol;
bc  = addBC([], f( top |  bot), 'flux'    , 0);
bc  = addBC(bc, f(~top & ~bot), 'pressure', 0);

xc  = .5*[1, 1, 1];
d   = sum(bsxfun(@minus, G.cells.centroids, xc).^2,2);
c   = find(d == min(d));
src = addSource([], c(1), 1*meter^3/day);

rock.perm = 100*milli*darcy*ones([G.cells.num,1]);
fluid     = initSingleFluid('mu',  100*centi*poise, ...
                            'rho', 1000*kilogram/meter^3);
state     = initState(G, [], 0);

S     = computeVirtualIP(G, rock, 2);
state = incompVEM(state, G, S, fluid, 'bc', bc, 'src', src);

plotCellData(G, state.pressure);