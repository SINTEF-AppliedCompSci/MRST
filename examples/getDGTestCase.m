function setup = getDGTestCase(name, varargin)

    setup = feval(lower(name), varargin);
    
end

%-------------------------------------------------------------------------%
function setup = simple1d(args) %#ok

    opt = struct('n', 100, 'nkr', 1, 'degree', [0,1,2]);
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
    tmodel = TransportModel(model, 'useTotalSaturation', true);
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
    setup = packSetup(state0, schedule, [], {modelFV}, {modelDG});
   
end

function setup = qfs_wo_2d(args) %#ok
    opt = struct('n', 50, 'nkr', 1, 'degree', [0,1,2]);
    [opt, discArgs] = merge_options(opt, args{:});

    G = computeGeometry(cartGrid([1,1]*opt.n, [500,500]*meter));
    G = createAugmentedGrid(G);
    G = computeCellDimensions(G);
    
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
                               'n'     , [2,2]                      );
    
    model  = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    pmodel = PressureModel(model);
    tmodel = TransportModel(model);
    tmodel.parentModel.useCNVConvergence = false;
    tmodel.parentModel.nonlinearTolerance = 1e-3;
    
    modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
    
    modelDG = cell(numel(opt.degree), 1);
    for dNo = 1:numel(opt.degree)
%         disc         = DGDiscretization(modelFV.transportModel, ...
%                                    'degree', opt.degree(dNo), discArgs{:});
        tmodelDG = TransportModelDG(model, 'formulation', 'missingPhase' , ...
                                           'degree'     , opt.degree(dNo), ...
                                           discArgs{:});
        tmodelDG.parentModel.useCNVConvergence = false;
        tmodelDG.parentModel.nonlinearTolerance = 1e-3;
        tmodelDG.parentModel.OutputStateFunctions = {};
        modelDG{dNo} = SequentialPressureTransportModel(pmodel, tmodelDG, 'parentModel', model);
    end

    time  = 2*year;
    rate  = sum(poreVolume(G, rock))/time;
    dt    = 20*day;
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
    
    setup = packSetup(state0, schedule, [], {{modelFV}}, {modelDG});
    
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