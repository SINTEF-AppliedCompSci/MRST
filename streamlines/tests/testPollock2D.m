%% Streamlines for a 2D Quater Five-Spot Well Pattern
% The quarter five-spot problem is one of the cases that are most widely
% used to test numerical methods. The setup consist of an injector and a
% producer located in diagonally oposite corners, with no-flow conditions
% along the model perimeter. This gives a symmetric flow pattern with high
% flow along the diagonal and stagnant points at the corners with no wells.

mrstModule add mimetic incomp streamlines

%% Construct model and compute flow field
G = cartGrid([25,25]);
G = computeGeometry(G);

rock.perm = repmat(10*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(0.3,            [G.cells.num, 1]);

fluid = initSimpleFluid('mu',  [1, 1]*centi*poise', ...
                        'rho', [1000, 1000]*kilogram/meter^3, ...
                        'n',   [2, 2]);

T = computeTrans(G, rock);
src = addSource([], [1, G.cells.num], [1.1, -1.1]);

x = initResSol(G, 0, 0);
x = incompTPFA(x, G, T, fluid, 'src', src);

%% Trace streamlines
% Pick start points on the diagonal orthogonal to the well-pair direction
% and trace streamlines forward and backward from these points. This will
% ensure maximal accuracy in the near-well regions where the streamlines
% are converging.
clf
h = plotGrid(G, 'facea', 0.3, 'edgea',0.1);
hold on;
cells = (G.cartDims(1):G.cartDims(1)-1:prod(G.cartDims)-1)';
streamline(pollock(G, x, cells));
x.flux = -x.flux;
streamline(pollock(G, x, cells, 'substeps', 1));
plot(G.cells.centroids(cells,1), G.cells.centroids(cells,2), ...
    'or','MarkerSize',8,'MarkerFaceColor',[.6 .6 .6]);
axis equal tight
hold off

%% Copyright notice

% #COPYRIGHT_EXAMPLE#