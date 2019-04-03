
verbose = false;
MODS = mrstModule;
mrstModule add mimetic mpfa incomp

%% Define and process geometry
% Construct a Cartesian grid of size 10-by-10-by-4 cells, where each cell
% has dimension 1-by-1-by-1. Because our flow solvers are applicable for
% general unstructured grids, the Cartesian grid is here represented using
% an unstructured formate in which cells, faces, nodes, etc. are given
% explicitly.
nx = 10; ny = 10;
G = cartGrid([nx, ny]);
G = twister(G, 0.1);
G = computeGeometry(G);
nc = G.cells.num;

%% Set rock and fluid data
% The only parameters in the single-phase pressure equation are the
% permeability $K$, which here is homogeneous, isotropic and equal 100 mD.
% The fluid has density 1000 kg/m^3 and viscosity 1 cP.
% We make a non diagonal rock tensor
rock = makeRock(G, 1e-3*darcy, 1);

% injcells  =  1 : 3;
% prodcells =  nc : -1 : (nc - 3);
injcells  =  1 : 3;
prodcells =  nc : -1 : (nc - 2);
radius    = 0.1;

rate = 1e-4*meter^3/day; 
bhp  = 1*barsa;

W = addWell([], G, rock, injcells, 'Comp_i', 1, 'Type', 'rate', 'Val', rate, ...
            'sign', 1, 'Radius', radius);
W = addWell(W, G, rock, prodcells, 'Comp_i', 1, 'Type', 'bhp', 'Val', bhp, ...
            'sign', -1, 'Radius', radius);

fluid = initSingleFluid('mu' , 1*centi*poise , ...
                        'rho', 1014*kilogram/meter^3);

gravity off
% 
% T_mpfa = computeMultiPointTrans(G, rock,'eta',1/3);
% resSol3 = initState(G, W, 0);

%% Solve MPFA pressure
% resSol3 = incompMPFA(resSol3, G, T_mpfa, fluid, 'wells', W);
% pressure = resSol3.pressure;

% clf
% plotCellData(G, pressure ./ barsa());
% title('Pressure: MPFA'); view(2), axis tight off
% colorbar('Location','SouthOutside');


% Compute the transmisibility matrix for mpfa
mpfastruct = computeMultiPointTrans2(G, rock, 'eta', 1/3);

%% Solve MPFA pressure
state = incompMPFA2(G, mpfastruct, 'wells', W);
pressure = state.pressure;

%% plotting

figure(1)
clf
plotCellData(G, pressure ./ barsa());
title('Pressure: MPFA'); view(2), axis tight off
colorbar('Location','SouthOutside');


IP = computeMimeticIP(G, rock, 'InnerProduct', 'ip_simple');
resSol2 = initState(G, W, 0);

%% Solve mimetic linear hybrid system
resSol2 = incompMimetic(resSol2, G, IP, fluid, 'wells', W);
pressure = resSol2.pressure;

figure(2)
clf
plotCellData(G, pressure ./ barsa());
title('Pressure: mimetic'); view(2), axis tight off
colorbar('Location','SouthOutside');