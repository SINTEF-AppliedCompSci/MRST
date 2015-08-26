function [state0, schedule, model] = setupSimpleOW()
    mrstModule add deckformat ad-fi ad-refactor ad-props

    fn = fullfile('SINTEF', 'simpleOW', 'simple10x1x10.data');

    [deck, schedule, model] = setupADcase(fn);

    gravity on

    state0 = initResSol(model.G, deck.PROPS.PVCDO(1), [.15, .85]);
end