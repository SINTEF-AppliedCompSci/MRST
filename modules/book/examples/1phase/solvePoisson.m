%% Solve the Poisson problem
% In the first example we use a small Cartesian grid
G = computeGeometry(cartGrid([5 5],[1 1]));

%% Define discrete operators
% Since we impose no-flow boundary conditions, we restrict the connections
% to the interior faces only
N = G.faces.neighbors;
N = N(all(N ~= 0, 2), :);
nf = size(N,1);
nc = G.cells.num;
C = sparse([(1:nf)'; (1:nf)'], N, ones(nf,1)*[-1 1], nf, nc);
grad = @(x) C*x;
div  = @(x) -C'*x;
figure; spy(C);

%% Set up and solve the problem
p = initVariablesADI(zeros(nc,1));
q = zeros(nc, 1);               % source term
q(1) = 1; q(nc) = -1;           % -> quarter five-spot

eq    = div(grad(p))+q;         % equation
eq(1) = eq(1) + p(1);           % make solution unique
p     = -eq.jac{1}\eq.val;      % solve equation
clf, plotCellData(G, p);
