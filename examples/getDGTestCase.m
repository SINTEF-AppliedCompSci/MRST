function setup = getDGTestCase(name, varargin)

    setup = feval(lower(name), varargin);
    
end

%-------------------------------------------------------------------------%
function setup = simple1d(args) %#ok

    opt = struct('n', 100, 'degree', [0,1,2]);
    [opt, discArgs] = merge_options(opt, args{:});

    G = computeGeometry(cartGrid([opt.n,1], [opt.n,1]*meter));
    G = createAugmentedGrid(G);
    G = computeCellDimensions(G);
    
    rock  = makeRock(G, 1, 1);
    fluid = initSimpleADIFluid('phases', 'WO' , ...
                               'rho'   , [1,1], ...
                               'mu'    , [1,1], ...
                               'n'     , [1,1]);
    
    model  = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    pmodel = PressureModel(model);
    tmodel = TransportModel(model);
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
%     bc    = [];
%     bc = fluxside(bc, G, 'left' ,  1, 'sat', [1,0]);
%     bc = fluxside(bc, G, 'right', -1, 'sat', [1,0]);
    
    schedule = simpleSchedule(dtvec, 'W', W);

    sW     = 0.0;
    state0 = initResSol(G, 1, [sW,1-sW]);
%     state0.flux = zeros(G.faces.num,1);
%     state0.flux(modelFV.parentModel.operators.internalConn) = 1;
%     state0.flux(modelFV.transportModel.operators.internalConn) = 1;
%     state0.flux(bc.face) = 1;
    
    setup = packSetup(state0, schedule, [], modelFV, {modelDG});
   
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