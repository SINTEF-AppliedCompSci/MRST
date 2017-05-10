mrstModule add ad-core ad-eor ad-blackoil ad-props blackoil-sequential matlab_bgl

G = computeGeometry(cartGrid([50,1,1]));
rock = makeRock(G, 1, 1);

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                              'rho'   , [1000, 800, 0.67]*kilogram/meter^3, ...
                              'phases', 'WOG', ...
                              'mu'    , [100, 500, 0.011]*centi*poise);

sres = 0.05;
fluid = addSolventProperties(fluid, 'mu', 0.014*centi*poise, 'rho', 1.8*kilogram/meter^3, ...
      'sOres_i', sres, 'sOres_m', 0.5*sres, 'sSGres_i', sres, 'sSGres_m', 0.5*sres, 'mixPar', 1);
      
fluid.krW = coreyPhaseRelpermAD(2, sres, 1, 4*sres);
fluid.krO = coreyPhaseRelpermAD(2, sres, 1, 4*sres);
fluid.krG = coreyPhaseRelpermAD(2, sres, 1, 4*sres);
fluid.krS = coreyPhaseRelpermAD(2, sres, 1, 4*sres);


T = 1*year;
pv = poreVolume(G, rock);
injRate = 1*sum(pv)/T;
W = [];

W = verticalWell(W, G, rock,           1, 1, [], 'comp_i', [0, 0, 0, 1], 'type', 'rate', 'val', injRate);
W = verticalWell(W, G, rock, G.cells.num, 1, [], 'comp_i', [0, 0, 0, 1], 'type', 'bhp' , 'val', 0      );

% fSt = G.faces.centroids(:,1) == min(G.faces.centroids(:,1));
% fEd = G.faces.centroids(:,1) == max(G.faces.centroids(:,1));
% bc = addBC([], find(fSt), 'flux', injRate, 'sat', [0 0 1]);
% bc = addBC(bc, find(fEd), 'pressure', 0, 'sat', [0 0 1]);

% model = ThreePhaseBlackOilModel(G, rock, fluid);


model = FourPhaseSolventModel(G, rock, fluid);

sres = 0.1;
state0 = initResSol(G, 50*barsa, [sres 1-3*sres sres, sres]);
state0.wellSol = initWellSolAD(W, model, state0);

nStep = 50;
dT = repmat(T/nStep, nStep, 1);
schedule = simpleSchedule(dT, 'W', W);

[ws, statesWOS, reports] = simulateScheduleAD(state0, model, schedule);


%%

W = [];

W = verticalWell(W, G, rock, 1, 1, [],           'comp_i', [1, 0], 'type', 'rate', 'val', injRate);
W = verticalWell(W, G, rock, G.cells.num, 1, [], 'comp_i', [1, 0], 'type', 'bhp', 'val', 0);


fluidWO = initSimpleADIFluid('n', [2, 2], ...
                             'rho', [1000, 800]*kilogram/meter^3, ...
                             'phases', 'WO', ...
                             'mu', [100, 500]*centi*poise);
                         

model = TwoPhaseOilWaterModel(G, rock, fluidWO);

state0 = initResSol(G, 0, [sres, 1-sres]);
state0.wellSol = initWellSolAD(W, model, state0);

nStep = 100;
dT = repmat(T/nStep, nStep, 1);
schedule = simpleSchedule(dT, 'W', W);

[ws, statesWO, reports] = simulateScheduleAD(state0, model, schedule);
