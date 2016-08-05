clc; clear all; close all

mrstModule add mimetic
mrstModule add mpfa
mrstModule add streamlines
addpath('../../vem/mat/VEM2D/')

xmax = 1;
ymax = 1;
nx   = 50;
ny   = 50;

G = cartGrid([nx,ny],[xmax,ymax]);

G.nodes.coords = twister(G.nodes.coords);
G = sortEdges(G);
G = computeVEM2DGeometry(G);

%%  Set BC
tol = 1e-6;
boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);
left = abs(G.faces.centroids(boundaryEdges,1)) < tol;
right= abs(G.faces.centroids(boundaryEdges,1) - xmax) < tol;
             
neuman = ~left& ~right;

bc_MRST = addBC([], boundaryEdges(left), 'pressure', 0);
bc_MRST = addBC(bc_MRST, boundaryEdges(right), 'pressure', 100);
bc_MRST = addBC(bc_MRST, boundaryEdges(neuman), 'flux', 0);

bc_VEM = VEM2D_addBC([], G, boundaryEdges(left), 'pressure', 0);
bc_VEM = VEM2D_addBC(bc_VEM, G, boundaryEdges(right), 'pressure', 100);
bc_VEM = VEM2D_addBC(bc_VEM, G, boundaryEdges(neuman), 'flux', 0);


%% Set fluid and rock properties
gravity reset off 

fluid = initSingleFluid('mu' , 1    , ...
                        'rho', 1);

rock.poro = ones(G.cells.num,1);

rock.perm = ones([G.cells.num,1]);



%% Initialize state
sInit = initState(G, [], 0, [0.0,1]);
S     = computeMimeticIP(G, rock, 'Verbose', true);
transTPFA = computeTrans(G,rock);
transMPFA = computeMultiPointTrans(G, rock);

%% Solve Laplace

%% TPFA
sTPFA = incompTPFA(sInit, G, transTPFA, fluid, 'bc',bc_MRST);

%% MIMETIC
sMIM  = incompMimetic(sInit, G, S, fluid,'bc',bc_MRST);

%% MPFA
sMPFA  = incompMPFA(sInit, G, transMPFA, fluid,'bc',bc_MRST);

%% VEM1
[sVEM1, G1] = VEM2D(G,0,bc_VEM,1,'cellAverages',true);

%% VEM2
sVEM2 = VEM2D(G,0,bc_VEM,2);
