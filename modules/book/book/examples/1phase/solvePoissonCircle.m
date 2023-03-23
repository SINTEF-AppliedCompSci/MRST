%% Grid 
G = cartGrid([20 20],[1 1]);
G = computeGeometry(G);
r1 = sum(bsxfun(@minus,G.cells.centroids,[0.5 1]).^2,2);
r2 = sum(bsxfun(@minus,G.cells.centroids,[0.5 0]).^2,2);
clf, plotCellData(G, double((r1>0.16) & (r2>0.16)) );

G = extractSubgrid(G, (r1>0.16) & (r2>0.16));

%% Grid information
N = G.faces.neighbors;
N = N(all(N ~= 0, 2), :);
nf = size(N,1);
nc = G.cells.num;

%% Operators
C = sparse([(1:nf)'; (1:nf)'], N, ...
    ones(nf,1)*[-1 1], nf, nc);
grad = @(x) C*x;
div  = @(x) -C'*x;
spy(C); set(gca,'XTick',[],'YTick',[]); xlabel([]);

%% Assemble and solve equations
p = initVariablesADI(zeros(nc,1));
q = zeros(nc, 1);             % source term
q(1) = 1; q(nc) = -1;         % -> quarter five-spot

eq    = div(grad(p))+q;       % equation
eq(1) = eq(1) + p(1);         % make solution unique
p     = -eq.jac{1}\eq.val;    % solve equation

%% Plot result
clf, plotCellData(G,p);
