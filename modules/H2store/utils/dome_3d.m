%% Load necessary MRST modules
mrstModule add deckformat

%% Set up basic grid parameters
dims = [10, 10, 10];
L = [50, 50, 10];
G = cartGrid(dims, L);
ub = 8;  % Upper boundary for deformation in z-direction

%% Apply deformation to the grid nodes in 3D
z = G.nodes.coords(:,3);
x = G.nodes.coords(:,1);
y = G.nodes.coords(:,2);
G.nodes.coords(:,3) = z - 1 * sin(pi * y / L(2)) .* sin(pi * x / L(1)) .* sin(pi * z / ub) .* (z < ub);
G = computeGeometry(G);

% Identify cells in the upper half of the grid along the y-axis
ind = G.cells.centroids(:,2) > (L(2) / 2);
% Optional visualization: plotGrid(G, ind, 'FaceColor', 'None'), view(3)

%% Define refined grid with padding and a custom depth
pad = 100;                  % Padding size for the grid extension
refine = 8;                 % Refinement factor
depth = -300;               % Custom depth for Z-coordinates
myres = [51, 51, 2];        % Grid resolution
Ldims = [5000, 5000, 100];  % Domain dimensions

%% Generate basic corner-point grid definition
grdecl = simpleGrdecl(myres, 0.0, 'undisturbed', true, 'flat', true);
grdecl.COORD(1:3:end) = grdecl.COORD(1:3:end) * Ldims(1);
grdecl.COORD(2:3:end) = grdecl.COORD(2:3:end) * Ldims(2);
grdecl.COORD(3:3:end) = grdecl.COORD(3:3:end) * Ldims(3);
grdecl.ZCORN = grdecl.ZCORN * Ldims(3) * 2;

%% Apply deformation using Gaussian bump in Z-coordinates
[X, Y, Z, lineIx] = buildCornerPtNodes(grdecl);
ZCORN = Z + 100 .* exp(-(X - 2500).^2 / (2 * 500^2) - (Y - 2500).^2 / (2 * 500^2));
grdecl.ZCORN = ZCORN + depth;
grdecl.ZCORN = reshape(grdecl.ZCORN, 1, []);

%% Apply padding and refinement to the grid
grdecl = padGrdecl(grdecl, [0 0 1], [0 0; 0 0; [pad pad]]);
grdecl = refineGrdecl(grdecl, [1 1 refine]);
G = processGRDECL(grdecl);
shift = 450;
G.nodes.coords(:,3) = G.nodes.coords(:,3) - shift;  % Adjust Z-coordinates by shifting downwards
G = computeGeometry(G);

%% Classify geological layers based on the Z-index
[i, j, k] = ind2sub(G.cartDims, G.cells.indexMap);
underrock = ismember(k, [1:8]);
bedrock = ismember(k, [9:12]);
toprock = ismember(k, [25:32]);
caprock = ismember(k, [20:24]);
aquifer = ~bedrock & ~caprock & ~underrock & ~toprock;

%% Visualize geological layers with color coding
yd = 2500;
plotGrid(G, caprock & G.cells.centroids(:,2) < yd, 'FaceColor', 'b');
plotGrid(G, bedrock & G.cells.centroids(:,2) < yd, 'FaceColor', 'r');
plotGrid(G, underrock & G.cells.centroids(:,2) < yd, 'FaceColor', 'g');
plotGrid(G, toprock & G.cells.centroids(:,2) < yd, 'FaceColor', 'm');
plotGrid(G, aquifer & G.cells.centroids(:,2) < yd);
title('Geological Layers of the Dome Aquifer');
view(0, 180);
