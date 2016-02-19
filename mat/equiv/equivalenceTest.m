addpath('../VEM2D/');

%   TEST 1: unit square

G = cartGrid([1,1], [2,2]);

X = G.nodes.coords;
X = bsxfun(@minus, X, [1,1]);
G.nodes.coords = X;

f = @(X) X(:,1);

k = 1;

G = computeVEM2DGeometry(G,f,k);
AK_FEM = FEM2D_loc(G,1, k);
AK_VEM = G.cells.AK{1};

err = norm(AK_FEM - AK_VEM, 'fro');

fprintf('Result from unit square comparison: \t  %d \n\n', err);

%   TEST 2: unit square

G = unitSquareTri([1,1],[2,2]);

X = G.nodes.coords;
X = bsxfun(@minus, X, [1,1]);
G.nodes.coords = X;

f = @(X) X(:,1);

k = 1;

G = computeVEM2DGeometry(G,f,k);

AK_FEM = FEM2D_loc(G,1,1);

AK_VEM = G.cells.AK{2};

err = norm(AK_FEM - AK_VEM, 'fro');

fprintf('Result from unit triangle comparison: \t  %d \n\n', err);

