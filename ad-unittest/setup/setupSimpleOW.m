function [state0, schedule, model] = setupSimpleOW()
    mrstModule add deckformat ad-fi ad-refactor

    fn = fullfile('SINTEF', 'simpleOW', 'simple10x1x10.data');

    [G, rock, fluid, deck, schedule] = setupADcase(fn);

    gravity on

    state0 = initResSol(G, deck.PROPS.PVCDO(1), [.15, .85]);
    model = TwoPhaseOilWaterModel(G, rock, fluid, 'inputdata', deck);
end