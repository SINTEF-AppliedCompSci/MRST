G = computeGeometry(cartGrid([100,1]));
rock = makeRock(G, repmat(100*milli*darcy,[G.cells.num,1]), ...
    repmat(0.2,[G.cells.num,1]));

fluid = initSimpleFluid('mu' , [   1,    1] .* centi*poise     , ...
                        'rho', [1000, 1000] .* kilogram/meter^3, ...
                        'n'  , [   2,    2]);
bc = fluxside([], G, 'Left',   1, 'sat', [1 0]);
bc = fluxside(bc, G, 'Right', -1, 'sat', [0 1]);

hT = computeTrans(G, rock);
rSol = initState(G, [], 0, [0 1]);

rSol = incompTPFA(rSol, G, hT, fluid, 'bc', bc);
rSole = explicitTransport(rSol, G, 10, rock, fluid, 'bc', bc);
rSoli = implicitTransport(rSol, G, 10, rock, fluid, 'bc', bc);
rSolt = rSol; n = 10;
for i=1:n
    rSolt = implicitTransport(rSolt, G, 10/n, rock, fluid, 'bc', bc);
end
plot(G.cells.centroids(:,1), rSole.s(:,1),'bo-', ...
    G.cells.centroids(:,1), rSoli.s(:,1), 'r*-', ...
    G.cells.centroids(:,1), rSolt.s(:,1), 'gs-');