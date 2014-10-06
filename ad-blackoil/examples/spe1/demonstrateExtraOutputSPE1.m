clear all
mrstModule add ad-fi deckformat mrst-gui ad-core ad-blackoil


[G, rock, fluid, deck, state] = setupSPE1();
model = selectModelFromDeck(G, rock, fluid, deck);


schedule = convertDeckScheduleToMRST(G, model, rock, deck);

schedule.step.val = schedule.step.val(1:5);
schedule.step.control = schedule.step.control(1:5);

model.extraStateOutput = true;
%% Run the entire schedule
[wellSols, states, report] = simulateScheduleAD(state, model, schedule);
%% Output surface fluxes from the stored properties
states = convertReservoirFluxesToSurface(model, states);