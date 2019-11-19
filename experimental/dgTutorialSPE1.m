mrstModule add dg ad-core ad-props ad-blackoil blackoil-sequential weno ...
    vem vemmech mrst-gui
mrstVerbose on

%%

[G, rock, fluid, deck, state0] = setupSPE1();
G = computeCellDimensions2(G);

model = selectModelFromDeck(G, rock, fluid, deck);

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(model, deck);

%%

modelSI = getSequentialModelFromFI(model);
% modelSI.transportModel.conserveWater = true;

%%

[wsFV, stFV, repFV] = simulateScheduleAD(state0, modelSI, schedule);

%%

transportModel = TransportBlackOilModelDG(G, rock, fluid, ...
                               'disgas', modelSI.transportModel.disgas, ...
                               'vapoil', modelSI.transportModel.vapoil, ...
                               'degree', 1);
modelDG        = SequentialPressureTransportModelDG(modelSI.pressureModel, transportModel);

[wsDG, stDG, repDG] = simulateScheduleAD(state0, modelDG, schedule);

%%

plotWellSols({wsDG, wsFV})

%%

sd = cellfun(@(s1, s2) compareStates(s1, s2), stFV, stDG, 'unif', false);
close all
plotToolbar(G, sd);