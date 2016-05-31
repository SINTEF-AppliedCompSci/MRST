clc; clear; close all;
%%  MAKE GRID

G = cartGrid([51,51],[1,1]);
G = computeGeometry(G);
        
%%  SET BCs

C = [0.5,0.5];
gD = @(X) -log(sqrt(sum(bsxfun(@minus,X,C).^2,2)))/(2*pi) + 10;

bcEdges = any(G.faces.neighbors == 0,2);
bc = addBC([], find(bcEdges), 'pressure', gD(G.faces.centroids(bcEdges ,:)));

%% FLUID AND ROCK PROPERTIES

gravity off 
fluid = initSingleFluid('mu' , 1, 'rho', 1);
rock  = struct('perm', ones(G.cells.num,1), 'poro', ones(G.cells.num,1));

%% ADD SOURCE

srcCell = find(all(abs(G.cells.centroids-C(1))<1e-6,2));
src = addSource([],srcCell,1);

%% INITIALIZE STATE

sInit1 = initState(G, [], 0);
sInit2 = initState(G, [], 0, [1,0]);
S      = computeMimeticIP(G, rock, 'Verbose', true);
T      = computeTrans(G,rock);

%% SOLVE

sol1 = incompTPFA(sInit1, G, T, fluid, 'bc', bc, 'src', src);
sol2 = incompTPFA(sInit2, G, T, fluid, 'bc', bc, 'src', src);

%% PLOT SOLUTIONS ALONG W.R.T. RADIUS

r = sqrt(sum(bsxfun(@minus,G.cells.centroids,C).^2,2));
xL = linspace(0.5,1,1000);
XL = [xL', xL'];
rL = sqrt(sum(bsxfun(@minus, XL, C).^2,2));

plot(rL,gD(XL))
set(gcf, 'DefaultTextInterpreter', 'Latex')
hold on
plot(r, sol1.pressure,'.', r, sol2.pressure,'.');

h = legend('Exact solution', 'sInit1', 'sInit2');
set(h, 'interpreter', 'latex');
xlabel('$r$'); ylabel('$p$')