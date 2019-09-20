% Test from Danilov Vassilevski 5.1 for testing the monotonicity.
% There is a monotonicity test in the Zhang Kobaisi paper, and the
% NTPFA fails.
%{

Summary
-------
Cart grid
k1=100,k2=10,k3=1
tpfa: monotone
ntpfa: not (strange because we don't do any correctHAPs)
mfd: not
However, all are monotone using k1=k2=k3=1

%}

clear all
close all

mrstModule add ad-core ad-props incomp mrst-gui postprocessing ...
    ntpfa-kobaisi mimetic

refine = 1;
n = 10;
N = [n,n,n]*refine;
G = cartGrid(N,[1,1,1]);
G = computeGeometry(G);

% remove [0.4,0.6]^3
ii = vecnorm(G.cells.centroids - [0.5,0.5,0.5],2,2) < 0.1;
%plotGrid(G, ii);
G = removeCells(G, ii);
G = computeGeometry(G);

% Rotation matrices
Rx = @(x) [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
Ry = @(x) [cos(x) 0 sin(x); 0 1 0; -sin(x) 0 cos(x)];
Rz = @(x) [cos(x) -sin(x) 0; sin(x) cos(x) 0; 0 0 1];

% Tensor
k1 = 100; 
k2 = 10;
k3 = 1;
thx = pi/3;
thy = pi/4;
thz = pi/6;
K = Rz(-thz)*Ry(-thy)*Rx(-thx)*diag([k1, k2, k3])*Rx(thx)*Ry(thy)*Rz(thz);
perm = [K(1,1), K(1,2), K(1,3), K(2,2), K(2,3), K(3,3)];
rock = makeRock(G, perm, 1);

fluid = initSingleFluid('mu', 1, 'rho', 1);

% bc
bf = boundaryFaces(G);
bc = [];
outer_faces = vecnorm(G.faces.centroids(bf,:) - [0.5,0.5,0.5],2,2) > 0.3;
outer_value = 2;
bc = addBC(bc, bf(outer_faces), 'pressure', outer_value);
inner_faces = ~outer_faces;
inner_value = 1;
bc = addBC(bc, bf(inner_faces), 'pressure', inner_value);
bcntpfa = convertBC2FlowNTPFA(G, bc);

titles = {};
states = {};

%% TPFA
state0 = initResSol(G, 0, 1);
T = computeTrans(G, rock);
state = incompTPFA(state0, G, T, fluid, 'bc', bc);
p0 = state.pressure;
titles{end+1} = 'tpfa';
states{end+1} = state;

%% ntpfa kobaisi
tol = 1e-14;
mrstVerbose on
interpFace=findHAP(G,rock,bcntpfa);
OSflux=findOSflux(G,rock,bcntpfa,interpFace);
state=FlowNTPFA(G,bcntpfa,fluid,[],OSflux,p0,tol,1000);
mrstVerbose off
titles{end+1} = 'ntpfa-kobaisi';
states{end+1} = state;

%% mfd
IP = computeMimeticIP(G,rock);
state = incompMimetic(state0,G,IP,fluid,'bc',bc);
titles{end+1} = 'mfd';
states{end+1} = state;

close all
for i = 1 : numel(states)
    disp(titles{i})
    disp([min(states{i}.pressure), max(states{i}.pressure)])
    figure
    plotToolbar(G, states{i}),view(3),colorbar
    title(titles{i})
end

