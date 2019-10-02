function setup = getDGTestCase(name, varargin)

    setup = feval(lower(name), varargin);
    
end

%-------------------------------------------------------------------------%
function setup = simple1d(args) %#ok

    opt = struct('n', 100, 'nkr', 1, 'degree', [0,1,2]');
    [opt, discArgs] = merge_options(opt, args{:});

    G = computeGeometry(cartGrid([opt.n,1], [opt.n,1]*meter));
    G = createAugmentedGrid(G);
    G = computeCellDimensions(G);
    [G.cells.equal, G.faces.equal] = deal(true);
    
    rock  = makeRock(G, 1, 1);
    fluid = initSimpleADIFluid('phases', 'WO' , ...
                               'rho'   , [1,1], ...
                               'mu'    , [1,1], ...
                               'n'     , [1,1].*opt.nkr);
    
    model  = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    pmodel = PressureModel(model);
    tmodel = TransportModel(model);
    tmodel.parentModel.useCNVConvergence = false;
    tmodel.parentModel.nonlinearTolerance = 1e-3;
    
    modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
    
    modelDG = cell(size(opt.degree,1), 1);
    for dNo = 1:size(opt.degree,1)
        disc         = DGDiscretization(modelFV, ...
                                   'degree', opt.degree(dNo,:), discArgs{:});
        tmodelDG     = TransportModelDG(model, 'disc', disc);
        tmodelDG.disc.limiter = getLimiter(tmodelDG, 'tvb', 0, 'plot', true);
        tmodelDG.parentModel.useCNVConvergence = false;
        tmodelDG.parentModel.nonlinearTolerance = 1e-3;
        tmodelDG.parentModel.OutputStateFunctions = {};
        modelDG{dNo} = SequentialPressureTransportModel(pmodel, tmodelDG, 'parentModel', model);
    end

    time  = 2*opt.n;
    dt    = opt.n/100;
    dtvec = rampupTimesteps(time, dt, 0);
    
    W = [];
    W = addWell(W, G, rock, 1    , 'type', 'bhp', 'val', opt.n-1, 'compi', [1,0], 'WI', 9999);
    W = addWell(W, G, rock, opt.n, 'type', 'bhp', 'val', 0, 'compi', [1,0], 'WI', 9999);
    
    schedule = simpleSchedule(dtvec, 'W', W);

    sW     = 0.0;
    state0 = initResSol(G, 1, [sW,1-sW]);
    setup = packSetup(state0, schedule, {{model}}, {{modelFV}}, {modelDG});
   
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function setup = qfs_wo_2d(args) %#ok
    opt = struct('n', 20, 'nkr', 1, 'degree', [0,1]');
    [opt, discArgs] = merge_options(opt, args{:});

    G = computeGeometry(cartGrid([1,1]*opt.n, [500,500]*meter));
    G = createAugmentedGrid(G);
%     G = computeCellDimensions(G);
    G = computeCellDimensions2(G);
    
    if 0
        rng(2019)
        perm = logNormLayers([opt.n, opt.n, 1])*10*milli*darcy;
        poro = perm./max(perm)*0.8;
    else
        perm = 100*milli*darcy;
        poro = 0.4;
    end
    rock  = makeRock(G, perm, poro);
    fluid = initSimpleADIFluid('phases', 'WO'                       , ...
                               'rho'   , [1000,800]*kilogram/meter^3, ...
                               'mu'    , [0.5,1]*centi*poise        , ...
                               'n'     , [1,1]*opt.nkr              );
    
    model  = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    pmodel = PressureModel(model);
    tmodel = TransportModel(model);
    tmodel.parentModel.useCNVConvergence = false;
    tmodel.parentModel.nonlinearTolerance = 1e-3;
    
    modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
    
    modelDG = cell(size(opt.degree,1), 1);
    for dNo = 1:numel(opt.degree)
        tmodelDG = TransportModelDG(model, 'formulation', 'totalSaturation' , ...
                                           'degree'     , opt.degree(dNo,:), ...
                                           discArgs{:});
        tmodelDG.parentModel.useCNVConvergence = false;
        tmodelDG.parentModel.nonlinearTolerance = 1e-3;
        tmodelDG.parentModel.OutputStateFunctions = {};
        modelDG{dNo} = SequentialPressureTransportModel(pmodel, tmodelDG, 'parentModel', model);
    end

    time  = 2*year;
    rate  = sum(poreVolume(G, rock))/time;
    dt    = 30*day;
    dtvec = rampupTimesteps(time, dt);
    
    W = [];
    W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate     , 'compi', [1,0]);
    W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 100*barsa, 'compi', [1,0]);
    
    schedule = simpleSchedule(dtvec, 'W', W);

    sW     = 0.0;
    state0 = initResSol(G, 100*barsa, [sW,1-sW]);
    setup = packSetup(state0, schedule, {{model}}, {{modelFV}}, {modelDG});
    
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function setup = qfs_wog_3d(args) %#ok
    
    gravity reset on

    opt = struct('n', 5, 'nkr', 2, 'degree', [0,1]', 'pebi', true);
    [opt, discArgs] = merge_options(opt, args{:});

    n = opt.n;
    l = 500*meter;
    d = 0;
    h = 100*meter;
    nl = 5;
    
    if opt.pebi
        rng(2019);
        G = pebiGrid(l/n, [1,1]*l, 'wellLines', {[d,d], [l-d, l-d]}, 'wellRefinement', false, 'wellGridFactor', 0.15);
    else
        G = cartGrid([n,n], [1,1]*l);
        G.cells.tag = false(G.cells.num,1);
        G.cells.tag([1;G.cells.num]) = true;
        ix = strcmpi(G.type, 'tensorGrid');
        G.type = G.type(~ix);
    end
    G = makeLayeredGrid(G, repmat(h/nl, nl, 1));
    G.cells.tag = repmat(G.cells.tag, nl, 1);
    
    G = computeGeometry(G);
    G = createAugmentedGrid(G);
    G = computeCellDimensions2(G);
    [G.cells.equal, G.faces.equal] = deal(false);
    
    perm = 100*milli*darcy;
    poro = 0.4;
    
    rock  = makeRock(G, perm, poro);
    fluid = initSimpleADIFluid('phases', 'WOG'                       , ...
                               'rho'   , [1000,800,500]*kilogram/meter^3, ...
                               'mu'    , [0.5,1,0.1]*centi*poise        , ...
                               'n'     , [1,1,1]*opt.nkr              );
    
    model  = GenericBlackOilModel(G, rock, fluid);
    pmodel = PressureModel(model);
    tmodel = TransportModel(model);
    tmodel.parentModel.useCNVConvergence = false;
    tmodel.parentModel.nonlinearTolerance = 1e-3;
    
    modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
    
    modelDG = cell(size(opt.degree,1), 1);
    for dNo = 1:size(opt.degree,1)
        tmodelDG = TransportModelDG(model, 'formulation', 'totalSaturation' , ...
                                           'degree'     , opt.degree(dNo,:), ...
                                           discArgs{:});
        tmodelDG.parentModel.useCNVConvergence = false;
        tmodelDG.parentModel.nonlinearTolerance = 1e-3;
        tmodelDG.parentModel.OutputStateFunctions = {};
        modelDG{dNo} = SequentialPressureTransportModel(pmodel, tmodelDG, 'parentModel', model);
    end

    time  = 2*year;
    rate  = 0.7*sum(poreVolume(G, rock))/time;
    dt    = 20*day;
    dtvec = rampupTimesteps(time, dt);
    
    W = [];
    W = addWell(W, G, rock, G.cells.tag & G.cells.centroids(:,1) < l/2, 'type', 'rate', 'val', rate, 'compi', [0,0,1]);
    W = addWell(W, G, rock, G.cells.tag & G.cells.centroids(:,1) > l/2, 'type', 'bhp' , 'val', 100*barsa, 'compi', [0,0,1]);
    
    schedule = simpleSchedule(dtvec, 'W', W);

    oil = G.cells.centroids(:,3) < h/3;
    state0 = initResSol(G, 100*barsa, [1,0,0]);
    state0.s(oil,:) = repmat([0,1,0], nnz(oil), 1);
    setup = packSetup(state0, schedule, {{model}}, {{modelFV}}, {modelDG});
   
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function setup = qfs_co2_small(args) %#ok
    % Model
    opt = struct('n', 20, 'degree', [0,1]);
    [opt, discArgs] = merge_options(opt, args{:});
    n = opt.n;
    G = cartGrid([n, n, 1], [1000, 1000, 10]);
    G = computeGeometry(G);
    G = createAugmentedGrid(G);
    G = computeCellDimensions2(G);
    K = 0.1*darcy;
    rock = makeRock(G, K, 0.25);

    fluid = initSimpleADIFluid('n', [2, 3],...
                               'phases', 'OG', ...
                               'rho', [100, 100]);
    [cf, info] = getCompositionalFluidCase('simple');

    model = GenericOverallCompositionModel(G, rock, fluid, cf, 'water', false);
    pmodel = PressureModel(model);
    tmodel = TransportModel(model);
    tmodel.parentModel.useCNVConvergence = false;
    tmodel.parentModel.nonlinearTolerance = 1e-3;
    
    modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
    
    modelDG = cell(numel(opt.degree), 1);
    for dNo = 1:numel(opt.degree)
        tmodelDG = TransportModelDG(model, 'formulation', 'totalSaturation' , ...
                                           'degree'     , opt.degree(dNo), ...
                                           discArgs{:});
        tmodelDG.parentModel.useCNVConvergence = false;
        tmodelDG.parentModel.nonlinearTolerance = 1e-3;
        tmodelDG.parentModel.OutputStateFunctions = {};
        modelDG{dNo} = SequentialPressureTransportModel(pmodel, tmodelDG, 'parentModel', model);
    end
    
    
    % Schedule
    time = 5*year;
    nstep = 100;

    irate = sum(model.operators.pv)/time;

    W = [];
    W = addWell(W, G, rock, 1, 'type', 'rate', 'val', irate, 'comp_i', [1, 0]);
    W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 25*barsa, 'comp_i', [1, 0]);

    for i = 1:numel(W)
        W(i).components = info.injection;
    end
    dt = rampupTimesteps(time, time/nstep);
    schedule = simpleSchedule(dt, 'W', W);
    %
    state0 = initCompositionalState(G, 50*barsa, info.temp, [1, 0], info.initial, model.EOSModel);
    
    setup = packSetup(state0, schedule, [], {{modelFV}}, {modelDG});
    
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function setup = spe1(args) %#ok

    opt = struct('degree', 0:2);
    [opt, discArgs] = merge_options(opt, args{:});
    
    [G, rock, fluid, deck, state0] = setupSPE1();
    G = computeCellDimensions2(G);

    model = selectModelFromDeck(G, rock, fluid, deck);

    pmodel = PressureModel(model);
    tmodel = TransportModel(model);
    tmodel.formulation = 'missingPhase';
    tmodel.parentModel.useCNVConvergence = false;
    tmodel.parentModel.nonlinearTolerance = 1e-3;
    modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);

    modelDG = cell(numel(opt.degree), 1);
    for dNo = 1:numel(opt.degree)
        disc         = DGDiscretization(modelFV, ...
                                   'degree', opt.degree(dNo), discArgs{:});
        tmodelDG     = TransportModelDG(model, 'disc', disc);
        tmodelDG.parentModel.useCNVConvergence = false;
        tmodelDG.parentModel.nonlinearTolerance = 1e-3;
        tmodelDG.parentModel.OutputStateFunctions = {};
        tmodelDG.formulation = 'missingPhase';
        modelDG{dNo} = SequentialPressureTransportModel(pmodel, tmodelDG, 'parentModel', model);
    end
    
    % Convert the deck schedule into a MRST schedule by parsing the wells
    schedule = convertDeckScheduleToMRST(model, deck);
    
    setup = packSetup(state0, schedule, {{model}}, {{modelFV}}, {modelDG});
    
end

%-------------------------------------------------------------------------%
function setup = packSetup(state0, schedule, modelFI, modelFV, modelDG, varargin)

    setup = struct('state0'      , state0  , ...
                   'schedule'    , schedule, ...
                   'modelFI'     , modelFI , ...
                   'modelFV'     , modelFV , ...
                   'modelDG'     , modelDG );      
    setup = merge_options(setup, varargin{:});
    
end