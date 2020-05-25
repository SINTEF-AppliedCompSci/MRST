function [description, state0, model, schedule, options, plotOptions] = validation_1d_compositional(varargin)
    description ...
        = ['Example from the compositional module with CO2 injection' ...
           'into 1D reservoir filled with Decane, CO2 and Methane'    ];
    if nargout == 1, return; end
    options = struct();
    pth = getDatasetPath('simplecomp');
    fn  = fullfile(pth, 'SIMPLE_COMP.DATA');
    % Read deck
    deck = readEclipseDeck(fn);
    deck = convertDeckUnits(deck);
    % Set up grid
    G = initEclipseGrid(deck);
    G = computeGeometry(G);
    % Set up rock
    rock  = initEclipseRock(deck);
    rock  = compressRock(rock, G.cells.indexMap);
    fluid = initDeckADIFluid(deck);
    % Define some surface densities
    fluid.rhoOS = 800;
    fluid.rhoGS = 10;
    % Define model
    eos = initDeckEOSModel(deck);
    model = GenericOverallCompositionModel(G, rock, fluid, eos.fluid, 'water', false);
    model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', false);
    % Get Schedule
    schedule = convertDeckScheduleToMRST(model, deck);
    % Manually set the injection composition
    [schedule.control.W.components] = deal([0, 1, 0]);
    % Injection is pure gas
    [schedule.control.W.compi] = deal([1, 0]);
    % The problem is defined at 150 degrees celsius with 75 bar initial
    % pressure. We set up the initial problem and make a call to the flash
    % routines to get correct initial composition.
    for i = 1:numel(schedule.control.W)
        schedule.control.W(i).lims = [];
    end
    % Initial conditions
    z0 = [0.6, 0.1, 0.3];
    T  = 150 + 273.15;
    p  = 75*barsa;
    state0 = initCompositionalState(G, p, T, 1, z0, eos);
    % Plotting
    plotOptions = {'plot1d'            , true          , ...
                   'lockCaxis'         , true          , ...
                   'Size'              , [800, 350]    , ...
                   'PlotBoxAspectRatio', [2.5,1,1]     , ...
                   'Projection'        , 'Orthographic', ...
                   'YLim'              , [0,1]         , ...
                   'View'              , [0,90]        };
end