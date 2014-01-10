%{
G = computeGeometry(cartGrid([30,30,1]));
src = [];
src = addSource(src, 1, 1);
src = addSource(src, G.cells.num, -1);
bc = [];
rock.perm = ones(G.cells.num, 1);
rock.poro = ones(G.cells.num, 1);

W = addWell([], G, rock, 1, 'type', 'bhp', 'val', 200*barsa);

x = initResSol(G, 0, 0);
x = mex_ifsh(x, G, rock, W, bc, src)

plotCellData(G, x.pressure, 'edgec', 'k', 'edgea', .1, 'facea', .625)
view(3), axis tight, grid on, colorbar southoutside
%}
%%{
clear
G = computeGeometry(cartGrid([3, 1]));
rock.perm = ones([G.cells.num, 1]);

%bc = pside([], G, 'left', 1);
bc = [];
bc = pside(bc, G, 'right', 0);
W  = [];
W  = addWell(W, G, rock, 1, 'type', 'bhp', 'val', 1);

x = initResSol(G, 0);

[x, wbhp, wflux] = mex_ifsh(x, G, rock, W, bc, []);
%}

x2 = initState(G, W, 0);

x2 = solveIncompFlow(x2, G, computeMimeticIP(G, rock),   ...
                     initSingleFluid('mu', 1, 'rho', 1), ...
                     'wells', W, 'bc', bc);

fprintf('Cell pressure error        : %12.5e [relative]\n', ...
        norm(x2.pressure - x.pressure, inf) / ...
        norm(x2.pressure             , inf));

fprintf('Face flux error            : %12.5e [relative]\n', ...
        norm(x2.flux - x.flux, inf) / norm(x2.flux, inf));

fprintf('Well perforation flux error: %12.5e [relative]\n', ...
        norm(vertcat(x2.wellSol.flux) - wflux, inf) / ...
        norm(vertcat(x2.wellSol.flux)        , inf));
