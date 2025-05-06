%% Grid 
load seamount
G = pebi(triangleGrid([x(:) y(:)], delaunay(x,y)));
G = computeGeometry(G);
clf, plotGrid(G);
axis tight off; set(gca,'XLim',[210.6 211.7]);

%% Grid information
N  = G.faces.neighbors;
N  = N(all(N ~= 0, 2), :);
nf = size(N,1);
nc = G.cells.num;

%% Operators
C = sparse([(1:nf)'; (1:nf)'], N, ...
    ones(nf,1)*[-1 1], nf, nc);
grad = @(x) C*x;
div  = @(x) -C'*x;
spy(C); set(gca,'XTick',[],'YTick',[]); xlabel([]);

%% Transmissibility
hT = computeTrans(G, struct('perm', ones(nc,1)));
cf = G.cells.faces(:,1);
T  = 1 ./ accumarray(cf, 1 ./ hT, [G.faces.num, 1]);
T  = T(all(N~=0,2),:);

%% Assemble and solve equations
p = initVariablesADI(zeros(nc,1));
q = zeros(nc, 1);
q([135 282 17]) = [-1 .5 .5];
eq    = div(T.*grad(p))+q;
eq(1) = eq(1) + p(1);
p     = -eq.jac{1}\eq.val;

%% Plot result
clf, plotCellData(G,p);
axis tight off; set(gca,'XLim',[210.6 211.7]);
