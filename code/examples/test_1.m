mrstModule add ad-core ad-eor ad-blackoil ad-props blackoil-sequential matlab_bgl

n = 10;
G = computeGeometry(cartGrid([n,n]));
rock = makeRock(G, 1, 1);

fluid = initSimpleADIFluid('n', [2, 2, 2], ...
                           'rho', [1, 1, 1], ...
                           'phases', 'WOG', ...
                           'mu', [1, 1, 1]*centi*poise);

fluid = addSolventProperties(fluid);

W = addWell([], G, rock, round(n*n/2), 'Comp_i', [0 0 0 1]);
state = initState(G, W, 0, [1,0,0,0]);
                       
model = OilWaterSolventModel(G, rock, fluid);
                       
