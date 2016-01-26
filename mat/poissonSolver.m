% mrstModule add mimetic;
require mimetic;

% Set up a Cartesian grid.
nx = 60;
ny = 60;
G = cartGrid([nx, ny], [1, 1]);
G = computeGeometry(G);
nc = G.cells.num;
nf = G.faces.num;
half_faces = G.cells.faces(:, 1);
nhf = numel(half_faces);

% Identify bottom and top faces, to apply the Dirichlet boundary conditions.
top_faces = (1 : G.faces.num)';
top_faces = top_faces(G.faces.centroids(:, 2) == 1);
top_p = ones(numel(top_faces), 1);

bottom_faces = (1 : G.faces.num)';
bottom_faces = bottom_faces(G.faces.centroids(:, 2) == 0);
bottom_p = 2*ones(numel(bottom_faces), 1);

dirich_faces = [top_faces; bottom_faces]; 
dirich_p = [top_p; bottom_p];

% Add source (roughly) in the middle
source_rhs = zeros(nc, 1);
source_cell = floor(nx/2 + nx*ny/2);
source_rhs(source_cell) = 0/G.cells.volumes(source_cell); 

% Deform the grid at top and bottom using the functions h and eta.
h = @(x) (ones(numel(x), 1) + 0.02*cos(6*pi*x));
eta = @(x) (0.1*cos(10*pi*x).*sin(pi*x));
x = G.nodes.coords;
xx = zeros(size(x));
x1 = x(:, 1);
xx(:, 1) = x1;
xx(:, 2) = -h(x1) + (eta(x1) + h(x1)).*x(:, 2);
G.nodes.coords = xx;

G = computeGeometry(G);

% We compute the mimetic scalar product
rock.perm = ones(G.cells.num, 1); % this is because MRST is a code for geosciences...
s   = computeMimeticIP(G, rock);
BI = s.BI; % Inverse of B (for B as defined on slide)

% We assemble the matrices C and D.
cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
C = - sparse(1:numel(cellNo), cellNo, 1);
D = sparse(1:numel(cellNo), double(half_faces), 1, numel(cellNo), G.faces.num);

% We identify the flux unknown corresponding to half faces where the Dirichlet
% condition is applied.
is_dirich_faces = false(nf, 1);
is_dirich_faces(dirich_faces) = true;

is_ext_faces  = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0) ;
is_int_faces = ~is_ext_faces;
nif = nnz(is_int_faces);
is_neumann_faces = is_ext_faces & ~is_dirich_faces;
is_neumann_half_faces = is_neumann_faces(half_faces);
neuman_half_faces = 1 : nhf;
neuman_half_faces = neuman_half_faces(is_neumann_half_faces);
nnhf = nnz(is_neumann_half_faces);

% We construct the right-hand side corresponding to the source term coming from the
% Dirichlet boundary conditons.
dirich_pii_rhs = zeros(nf, 1); % called pii instead of pi
dirich_pii_rhs(dirich_faces) = dirich_p;
dirich_rhs = -D*dirich_pii_rhs; 

% Reduce the system to the unknown variables.
D = D(:, is_int_faces);
N = - sparse(neuman_half_faces, 1 : nnhf , 1, nhf, nnhf);

% Assemble the source term and compute the right-hand side.
source_rhs = [zeros(nhf, 1); source_rhs; zeros(nif, 1); zeros(nnhf, 1)];

dirich_rhs = [dirich_rhs; zeros(nc + nif + nnhf, 1)];
rhs = source_rhs + dirich_rhs;

% Schur reduction

R = [[-C'; -D'; -N']*BI, eye(nc + nif + nnhf)];
A = [[C, D, N]; zeros(nc + nif + nnhf)];
Q = R*A; 
rhs = R*rhs;

% Solve the system.
sol = Q\rhs;

% Recover the potential.
p = sol([true(nc, 1); false(nif, 1); false(nnhf, 1)]);

% Plot the potential.
figure(2); clf;
plotCellData(G, p);
colorbar;

