
%   TEST 1: unit square

addpath('./..')

G = cartGrid([1,1], [2,2]);

X = G.nodes.coords;
X = bsxfun(@minus, X, [1,1]);
G.nodes.coords = X;
k = 1;

f = @(X) zeros(size(X,1),1);

G = computeVEM2DGeometry(G,f,k,6);

AK_FEM = [2/3, -1/6, -1/3, -1/6 ; ...
         -1/6, 2/3, -1/6, -1/3  ;
         -1/3, -1/6, 2/3, -1/6  ;
         -1/6, -1/3, -1/6, 2/3];

% load('FEM2D_1st_sq');
AK_VEM = G.cells.AK{1};

err = norm(AK_FEM - AK_VEM, 'fro');

fprintf('Result from unit square comparison: \t  %d \n\n', err);