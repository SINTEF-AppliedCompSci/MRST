%% Streamlines for a 3D Quater Five-Spot Well Pattern
% The setup consist of a vertical injector and a vertical producer located
% in the southwest and northeast corners, respectively.  This gives a
% symmetric areal flow pattern with high flow along the diagonal and
% stagnant points at the southeast and northwest corners. The setup is
% almost the same as in the <matlab:edit('testPollock2D.m') 2D test>,
% except that we now use a consistent discretization from the |mimetic|
% module.

mrstModule add mimetic incomp streamlines

%% Set up model
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

%% Trace streamline
clf
h = plotGrid(G, 'facea', 0.3, 'edgea',0.1);
hold on;
cells = (G.cartDims(1):G.cartDims(1)-1:prod(G.cartDims)/2-1)';
cells = [cells; cells+prod(G.cartDims)/2];
streamline(pollock(G, x, cells));
x.flux = -x.flux;
streamline(pollock(G, x, cells, 'substeps', 1));
plot3(G.cells.centroids(cells,1), G.cells.centroids(cells,2), ...
     G.cells.centroids(cells,3), ...
    'or','MarkerSize',8,'MarkerFaceColor',[.6 .6 .6]);
hold off
axis tight, view(30,40)

%% Copyright notice

% #COPYRIGHT_EXAMPLE#