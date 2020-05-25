function  [description, state0, model, schedule, options, plotOptions] = norne_field_bo(varargin)
    % Warning: This example will likely not run without tweaking the
    % solvers (see 'fieldModelNorneExample').
    description ...
        = ['The full Norne field model. Made available through the ', ...
           'Open Porous Media Project, https://opm-project.org/'    ];
    if nargout == 1, return; end

    require ad-core ad-blackoil ad-props deckformat
    options = struct('G_viz', []); % Sneak in grid for visualization in options
    
    gravity reset on
    mrstVerbose on
    useMex = true;
    opm = mrstPath('opm-tests');
    assert(~isempty(opm), 'You must register https://github.com/opm/opm-tests as a module!');
    [deck, output] = getDeckOPMData('norne', 'NORNE_ATW2013');
    % Build model from EGRID/INIT files
    egrid = readEclipseOutputFileUnFmt([output.opm.location, '.EGRID']);
    init = readEclipseOutputFileUnFmt([output.opm.location, '.INIT']);
    [Ge, rock_ecl, Ne, Te] = initGridFromEclipseOutput(init, egrid, 'outputSimGrid', true);
    G_viz = computeGeometry(Ge{1});
    G_sim = Ge{2};
    options.G_viz = G_viz;


    rock = initEclipseRock(deck);
    rock = compressRock(rock, G_sim.cells.indexMap);

    fluid = initDeckADIFluid(deck, 'useMex', useMex);
    % Setup model, but skip setting up the operators since we do not have a
    % proper grid
    model = GenericBlackOilModel(G_sim, [], fluid, 'disgas', true, 'vapoil', true, 'inputdata', deck);
    % Finally set up the connectivity graph from output
    model.rock = rock;
    model.operators = setupOperatorsTPFA(G_sim, rock, 'deck', deck, 'neighbors', Ne, 'trans', Te);
    % Set up everything
    [state0, model, schedule] = initEclipseProblemAD(deck, 'model', model, ...
                                                            'TimestepStrategy', ...
                                                            'ds', 'useCPR', true, ...
                                                            'useMex', useMex);

    model.toleranceCNV = 1e-2;
    model.toleranceMB = 1e-7;

    % Set well tolerances
    model.FacilityModel = ExtendedFacilityModel(model);
    model.FacilityModel.toleranceWellBHP = 1e-3;
    model.FacilityModel.toleranceWellRate = 5e-3;

    % Reset just in case
    model.FluxDiscretization = [];
    model.FlowPropertyFunctions = [];
    model.PVTPropertyFunctions = [];

    model.AutoDiffBackend.useMex = true;
    model.AutoDiffBackend.rowMajor = true;

    model = model.validateModel();
    % Use the alternative more rigorous crossflow definition for component
    % fluxes
    xflow = WellComponentTotalVolumeBalanceCrossflow(model);
    xflow.onlyLocalDerivatives = true;

    model.FacilityModel.FacilityFluxDiscretization.ComponentTotalFlux = xflow;

    useFlag = false;
    model.PVTPropertyFunctions.Viscosity.useSaturatedFlag = useFlag;
    model.PVTPropertyFunctions.ShrinkageFactors.useSaturatedFlag = useFlag;

    RvMax = model.PVTPropertyFunctions.getStateFunction('RvMax');
    RvMax.rvReduction = 0.5; % Not tested, but in the deck
    model.PVTPropertyFunctions = model.PVTPropertyFunctions.setStateFunction('RvMax', RvMax);
    
    plotOptions = {};

end
