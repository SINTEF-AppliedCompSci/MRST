%% Basic Grid Operations and Manipulations
% In this tutorial, we describe how you can use MRST to construct Cartesian
% and rectilinear grids for rectangular and non-rectangular domains and
% show how you can populate your grid with petrophysical properties. To
% this end, we will grid a rectangular domain, from which we cut out a half
% circle. This example is discussed in more detail in 
% <http://www.sintef.no/projectweb/mrst/jolts/ Jolt 2>

%% Tensor grid
clf
dx = 1 - 0.5*cos((-1:0.1:1)*pi);
x = -1.15 + 0.1*cumsum(dx);
y = 0:0.05:1;
G = tensorGrid(x, y.^0.75);
plotGrid(G, 'FaceColor', 'none');


%% Compute geometry
disp('G.cells = '); disp(G.cells)
G = computeGeometry(G);
disp('G.cells = '); disp(G.cells)
cla, plotCellData(G, G.cells.volumes);


%% Define indicator for circular sector
r = sqrt(sum(G.cells.centroids .^ 2, 2));
cla, plotCellData(G, double(r > 0.6));


%% Extract subgrid
disp('G.cells = '); disp(G.cells)
Gs = extractSubgrid(G, r > 0.6);

disp('Gs.cells = '); disp(Gs.cells)
cla, plotGrid(Gs, 'FaceColor', 'none');


%% Generate permeability
K = convertFrom(logNormLayers([G.cartDims, 1], 1), milli*darcy);
rock.perm = reshape(K(Gs.cells.indexMap), [], 1);
cla, plotCellData(Gs, rock.perm);


%% Perturbe inner nodes
fI = false([Gs.faces.num, 1]);
fI(boundaryFaces(Gs)) = true;

nI = true([Gs.nodes.num, 1]);
n  = mcolon(Gs.faces.nodePos(fI), Gs.faces.nodePos(fI) + 1);
nI(Gs.faces.nodes(n)) = false;

Gs.nodes.coords(nI,:) = Gs.nodes.coords(nI,:) + 0.007*randn([sum(nI), 2]);
cla, plotCellData(Gs, rock.perm);

%% Copyright Notice
%
% #COPYRIGHT_EXAMPLE#