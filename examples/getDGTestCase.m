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
    
    modelFV = deal(TransportOilWaterModel(G, rock, fluid));
    modelDG = cell(numel(opt.degree), 1);
    
    for dNo = 1:numel(opt.degree)
        disc         = DGDiscretization(modelFV, ...
                                   'degree', opt.degree(dNo), discArgs{:});
        modelDG{dNo} = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);
    end

    time  = opt.n;
    dt    = opt.n/100;
    dtvec = rampupTimesteps(time, dt, 0);
    bc    = [];
    bc = fluxside(bc, G, 'left' ,  1, 'sat', [1,0]);
    bc = fluxside(bc, G, 'right', -1, 'sat', [1,0]);
    
    schedule = simpleSchedule(dtvec, 'bc', bc);

    sW     = 0.0;
    state0 = initResSol(G, 100*barsa, [sW,1-sW]);
    state0.flux = zeros(G.faces.num,1);
    state0.flux(modelFV.operators.internalConn) = 1;
    state0.flux(bc.face) = 1;
    
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