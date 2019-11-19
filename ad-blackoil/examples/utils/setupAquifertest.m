function  [model, initState, schedule] = setupAquifertest(fn)

    deck = readEclipseDeck(fn);
    deck = convertDeckUnits(deck);

    % G = initEclipseGrid(deck, 'SplitDisconnected', false);
    G = initEclipseGrid(deck);
    G = computeGeometry(G);

    rock  = initEclipseRock(deck);
    rock  = compressRock(rock, G.cells.indexMap);
    % rock.perm = 1e8*rock.perm;

    % end-point scaling, so include grid
    fluid = initDeckADIFluid(deck, 'G', G);

    output = processAquifer(deck, G);
    aquifers     = output.aquifers;
    aquind       = output.aquind; 
    initval      = output.initval;
    aquiferprops = output.aquiferprops;

    model = AquiferBlackOilModel(G, rock, fluid, aquifers, aquind, aquiferprops, ...
                                 'modeltype', 'oilwater');

    if isfield(deck.RUNSPEC, 'VAPOIL')
        model.vapoil = deck.RUNSPEC.VAPOIL;
        model.vapoil = deck.RUNSPEC.VAPOIL;
    else
        model.vapoil = 0;        
        model.vapoil = 0;        
    end
    if isfield(deck.RUNSPEC, 'DISGAS')
        model.disgas = deck.RUNSPEC.DISGAS;
        model.disgas = deck.RUNSPEC.DISGAS;
    else
        model.disgas = 0;        
        model.disgas = 0;        
    end
    model.FacilityModel = selectFacilityFromDeck(deck, model);
    % nc = G.cells.num;
    % initState.pressure = 250*barsa*ones(nc, 1);
    % s = 0.3*ones(nc, 2);
    % s(:, 2) = 1 - s(:, 1);
    % initState.s = s;
    initState = initStateDeck(model, deck);
    initState.aquiferpressures = initval.pressures;
    initState.aquifervolumes   = initval.volumes;

    schedule = convertDeckScheduleToMRST(model, deck);

end