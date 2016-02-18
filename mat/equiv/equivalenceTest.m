G = cartGrid([1,1], [2,2]);

X = G.nodes.coords;
X = bsxfun(@minus, X, [1,1]);
G.nodes.coords = X;
G = computeGeometry(G);

G = mrstGridWithFullMappings(G);

AK_FEM = FEM2D_loc(G,1);

f = @(X) X(:,1);
[AK_VEM, ~, ~, ~] = VEM2D_loc(G,1,f);