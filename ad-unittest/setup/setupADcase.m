function [deck, schedule, model, rock] = setupADcase(fn)
    moddir = mrstPath('query', 'ad-testdata');
    fn = fullfile(moddir, fn);
    if ~exist(fn, 'file')
        error(['Did not find dataset at expected location: (', fn , ')'])
    end

    deck = readEclipseDeck(fn);
    deck = convertDeckUnits(deck);

    G = initEclipseGrid(deck);
    G = computeGeometry(G);

    rock  = initEclipseRock(deck);
    rock  = compressRock(rock, G.cells.indexMap);

    fluid = initDeckADIFluid(deck);
    
    model = selectModelFromDeck(G, rock, fluid, deck);
    schedule = convertDeckScheduleToMRST(G, model, rock, deck);
    
end