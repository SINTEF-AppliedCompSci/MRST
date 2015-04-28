try
    require mimetic incomp streamlines
catch
    mrstModule add mimetic incomp streamlines
end
G = cartGrid([25,25,2]);
G = computeGeometry(G);
rock.perm = repmat(10*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(0.3,            [G.cells.num, 1]);

fluid = initSimpleFluid('mu',  [1, 1]*centi*poise', ...
                        'rho', [1000, 1000]*kilogram/meter^3, ...
                        'n',   [2, 2]);

IP = computeMimeticIP(G, rock);
src = addSource([], [1, G.cells.num], [1, -1]);
x = initResSol(G, 0, 0);
x = incompMimetic(x, G, IP, fluid, 'src', src);

h = plotGrid(G, 'facea', 0.3, 'edgea',0.1);
hold on;
cells = (G.cartDims(1):G.cartDims(1)-1:prod(G.cartDims)-1)';
streamline(pollock(G, x, cells));
x.flux = -x.flux;
streamline(pollock(G, x, cells, 'substeps', 1));

hold off
axis equal tight
