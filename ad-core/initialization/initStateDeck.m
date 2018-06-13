function state = initStateDeck(model, deck)
    regions = getInitializationRegionsDeck(model, deck);
    state = initStateBlackOilAD(model, regions);
end