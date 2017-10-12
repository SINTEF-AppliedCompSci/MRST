function fluid = addPolymerProperties(fluid, varargin)
    polyDeck  = readEclipseDeck('POLYMER.DATA');
    polyDeck  = convertDeckUnits(polyDeck);
    polyFluid = initDeckADIFluid(polyDeck);

    fns = {'muWMult','dps','rrf','rhoR','adsInx','adsMax','ads',...
        'mixPar','cmax'};
    for i=1:numel(fns)
        fluid.(fns{i}) = polyFluid.(fns{i});
    end
    fluid.dps = 0;
    fluid.rhoR = 0;
end