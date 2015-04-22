function [G, rock, fluid, deck, state] = setupSPE9()
    % Read and process file.
    pth = getDatasetPath('spe9');
    fn  = fullfile(pth, 'BENCH_SPE9.DATA');

    deck = readEclipseDeck(fn);

    % The deck is given in field units, MRST uses metric.
    deck = convertDeckUnits(deck);

    G = initEclipseGrid(deck);
    G = computeGeometry(G);

    rock  = initEclipseRock(deck);
    rock  = compressRock(rock, G.cells.indexMap);

    % Create a special ADI fluid which can produce differentiated fluid
    % properties.
    fluid = initDeckADIFluid(deck);

    % The case includes gravity
    gravity reset on

    p0  = deck.SOLUTION.PRESSURE;
    sw0 = deck.SOLUTION.SWAT;
    sg0 = deck.SOLUTION.SGAS;
    s0  = [sw0, 1-sw0-sg0, sg0];
    % Gas in oil phase
    if isfield(deck.SOLUTION, 'RS')
        rs0 = deck.SOLUTION.RS;
    else
        rs0 = 0;
    end
    % Oil in gas phase
    if isfield(deck.SOLUTION, 'RV')
        rv0 = deck.SOLUTION.RV;
    else
        rv0 = 0;
    end
    
    state = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);
end