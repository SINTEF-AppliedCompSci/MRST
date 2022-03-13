function setup = getDGTestCase(name, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    setup = feval(lower(name), varargin);
    
end

%-------------------------------------------------------------------------%
function setup = simple1d(args) %#ok

    opt = struct('n', 100, 'nkr', 1, 'degree', {{[0,0], [1,0], [2,0], [3,0]}});
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
    fluid = restrictRelperms(fluid);
    
    model  = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    pmodel = PressureModel(model);
    tmodel = TransportModel(model);
    tmodel.parentModel.useCNVConvergence = false;
    tmodel.parentModel.nonlinearTolerance = 1e-3;
    
    modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
    
    modelDG = cell(numel(opt.degree), 1);
    for dNo = 1:numel(opt.degree)
        disc         = DGDiscretization(G, ...
                                   'degree'  , opt.degree{dNo}, discArgs{:});
        tmodelDG     = TransportModelDG(model, 'discretization', disc);
        tmodelDG.parentModel.useCNVConvergence = false;
        tmodelDG.parentModel.nonlinearTolerance = 1e-3;
        tmodelDG.storeUnlimited = true;
        modelDG{dNo} = SequentialPressureTransportModel(pmodel, tmodelDG, 'parentModel', model);
    end

    time  = 2*opt.n;
    dt    = opt.n/50;
    dtvec = rampupTimesteps(time, dt, 0);
    
    W = [];
    W = addWell(W, G, rock, 1    , 'type', 'bhp', 'val', opt.n-1, 'compi', [1,0], 'WI', 9999);
    W = addWell(W, G, rock, opt.n, 'type', 'bhp', 'val', 0, 'compi', [1,0], 'WI', 9999);
    
    schedule = simpleSchedule(dtvec, 'W', W);

    sW     = 0.0;
    state0 = initResSol(G, 1, [sW,1-sW]);
    setup = packSetup(state0, schedule, {{model}}, {{modelFV}}, {modelDG});
    setup.plot1d = true;
    
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function setup = qfs_wo_2d(args) %#ok
    opt = struct('n', 20, 'nkr', 1, 'degree', {{0,1,2}}, 'k', {{[]}}, 'useGenericFV', true, 'rotate', false);
    [opt, discArgs] = merge_options(opt, args{:});

    if isempty(opt.k{1})
        k = cell(numel(opt.degree),1);
        [k{:}] = deal([]);
        opt.k = k;
    else
        opt.degree = cellfun(@(k) max(sum(k,2)), opt.k, 'UniformOutput', false);
    end

    G = cartGrid([1,1]*opt.n, [500,500]*meter);
    if opt.rotate
        R = @(t) [cos(t), -sin(t); sin(t), cos(t)];
        G.nodes.coords = (R(pi/4)*(G.nodes.coords - 500/2)' + 500/2)';
    end
    G = computeGeometry(G);
    G = createAugmentedGrid(G);
    G = computeCellDimensions(G);
    G.cells.equal = true;
    
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
                               'mu'    , [1,1]*centi*poise        , ...
                               'n'     , [1,1]*opt.nkr              );
    fluid = restrictRelperms(fluid);
    
    mg = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    if opt.useGenericFV
        model  = mg;
        pmodel = PressureModel(model);
        tmodel = TransportModel(model, 'formulation', 'missingPhase');
        tmodel = TransportModel(model);
        tmodel.parentModel.useCNVConvergence = false;
        tmodel.parentModel.nonlinearTolerance = 1e-3;
    else
        model  = TwoPhaseOilWaterModel(G, rock, fluid);
        pmodel = PressureOilWaterModel(G, rock, fluid);
        tmodel = TransportOilWaterModel(G, rock, fluid);
        tmodel.useCNVConvergence = false;
        tmodel.nonlinearTolerance = 1e-3;
    end
    
    modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
    
    modelDG = cell(numel(opt.degree), 1);
    for dNo = 1:numel(opt.degree)
        tmodelDG = TransportModelDG(mg, 'formulation', 'totalSaturation', ...
                                        'degree'     , opt.degree{dNo}  , ...
                                        'k'          , opt.k{dNo}       , ...
                                        discArgs{:}                     );
        tmodelDG.parentModel.useCNVConvergence = false;
        tmodelDG.parentModel.nonlinearTolerance = 1e-3;
        tmodelDG.parentModel.OutputStateFunctions = {};
        modelDG{dNo} = SequentialPressureTransportModel(pmodel, tmodelDG, 'parentModel', model);
    end

    time  = 2*year;
    rate  = 1.5*sum(poreVolume(G, rock))/time;
    dt    = 10*day;
%     dtvec = rampupTimesteps(time, dt);
    dtvec = rampupTimesteps(time, dt, 0);
    
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
function setup = viscous_fingers(args) %#ok
    opt = struct('n', 50, 'nkr', 1, 'degree', {{0,1,2}}, 'k', {{[]}}, 'rotate', false);
    [opt, discArgs] = merge_options(opt, args{:});

    if isempty(opt.k{1})
        k = cell(numel(opt.degree),1);
        [k{:}] = deal([]);
        opt.k = k;
    else
        opt.degree = cellfun(@(k) max(sum(k,2)), opt.k, 'UniformOutput', false);
    end

    G = cartGrid([1,1/5]*opt.n, [500,100]*meter);
    if opt.rotate
        R = @(t) [cos(t), -sin(t); sin(t), cos(t)];
        G.nodes.coords = (R(pi/4)*(G.nodes.coords - 500/2)' + 500/2)';
    end
    G = computeGeometry(G);
    G = createAugmentedGrid(G);
    G = computeCellDimensions(G);
    G.cells.equal = true;
    
    if 0
        rng(2019)
        perm = logNormLayers(G.cartDims)*10*milli*darcy;
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
    fluid = restrictRelperms(fluid);
    
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    pmodel = PressureModel(model);
    pmodel = PressureOilWaterModel(G, rock, fluid);
    tmodel = TransportModel(model);
    tmodel.parentModel.useCNVConvergence = false;
    tmodel.parentModel.nonlinearTolerance = 1e-3;
    
    modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
    
    modelDG = cell(numel(opt.degree), 1);
    for dNo = 1:numel(opt.degree)
        tmodelDG = TransportModelDG(model, 'formulation', 'totalSaturation', ...
                                        'degree'     , opt.degree{dNo}  , ...
                                        'k'          , opt.k{dNo}       , ...
                                        discArgs{:}                     );
        tmodelDG.parentModel.useCNVConvergence = false;
        tmodelDG.parentModel.nonlinearTolerance = 1e-3;
        tmodelDG.parentModel.OutputStateFunctions = {};
        modelDG{dNo} = SequentialPressureTransportModel(pmodel, tmodelDG, 'parentModel', model);
    end

    time  = 2*year;
    rate  = sum(poreVolume(G, rock))/time;
    dt    = 10*day;
    dtvec = rampupTimesteps(time, dt, 0);
    
    
    bc = [];
    
    if 1
        bc = fluxside(bc, G, 'left', rate, 'sat', [1,0]);
        bc = pside(bc, G, 'right', 10*barsa, 'sat', [0,1]);
    else
        bc = fluxside(bc, G, 'left', rate, 'sat', [1,0]);
        f = bc.face(G.faces.centroids(bc.face,2) > -68);
        bc = addBC([], f, 'flux', rate/2, 'sat', [1,0]);
        bcp = pside([], G, 'right', 10*barsa, 'sat', [0,1]);
        f = bcp.face(G.faces.centroids(bcp.face,2) > 285);
        bc = addBC(bc, f, 'pressure', 10*barsa, 'sat', [1,0]);
    end
    schedule = simpleSchedule(dtvec, 'bc', bc);

    sW     = 0.0;
    state0 = initResSol(G, 100*barsa, [sW,1-sW]);
    setup = packSetup(state0, schedule, {{model}}, {{modelFV}}, {modelDG});
    
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function setup = qfs_wog_3d(args) %#ok
    
    gravity reset on

    opt = struct('n', 7, 'nkr', 2, 'degree', [0,1]', 'pebi', true);
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
    G = computeCellDimensions(G);
    [G.cells.equal, G.faces.equal] = deal(false);
    
    perm = 100*milli*darcy;
    poro = 0.4;
    
    rock  = makeRock(G, perm, poro);
    fluid = initSimpleADIFluid('phases', 'WOG'                          , ...
                               'rho'   , [1000,800,500]*kilogram/meter^3, ...
                               'mu'    , [0.5,1,0.1]*centi*poise        , ...
                               'n'     , [1,1,1]*opt.nkr              );
    fluid = restrictRelperms(fluid);
    
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
function setup = spe10_wo(args) %#ok
    
    opt = struct('layers', 10, 'dt', 15*day, 'T', 2000*day, 'I', 1:60, 'J', 1:220, 'degree', (0:1)');
    [opt, discArgs] = merge_options(opt, args{:});
    
    make2d = numel(opt.layers) == 1;
    [state0, model, schedule] = setupSPE10_AD('layers', opt.layers, 'I', opt.I, 'J', opt.J, 'make2d', make2d);

    fluid = restrictRelperms(model.fluid);
    
    G = model.G;
    G = computeGeometry(G);
    G = createAugmentedGrid(G);
    G = computeCellDimensions2(G);
    [G.cells.equal, G.faces.equal] = deal(true);
    
    model  = GenericBlackOilModel(G, model.rock, fluid, 'gas', false);
    pmodel = PressureModel(model);
    tmodel = TransportModel(model);
    tmodel.formulation = 'totalSaturation';
    tmodel.parentModel.useCNVConvergence = false;
    tmodel.parentModel.nonlinearTolerance = 1e-3;
    
    nls = NonLinearSolver('useLineSearch', true);
    
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

    setup = packSetup(state0, schedule, {{model}}, {{modelFV}}, {modelDG});
   
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function setup = qfs_co2_2d(args) %#ok
    % Model
    opt = struct('n', 20, 'degree', [0,1], 'useOverall', true, 'useGeneric', true);
    [opt, discArgs] = merge_options(opt, args{:});
    n = opt.n;
    G = cartGrid([n, n, 1], [1000, 1000, 10]);
    G = computeGeometry(G);
    G = createAugmentedGrid(G);
    G = computeCellDimensions(G);
    K = 0.1*darcy;
    rock = makeRock(G, K, 0.25);

    fluid = initSimpleADIFluid('n'     , [2, 3]    ,...
                               'phases', 'OG'      , ...
                               'rho'   , [100, 100]);
    fluid = restrictRelperms(fluid);
    [cf, info] = getCompositionalFluidCase('simple');

    if opt.useOverall
        if opt.useGeneric
            model  = GenericOverallCompositionModel(G, rock, fluid, cf, 'water', false);
            pmodel = PressureOverallCompositionModel(G, rock, fluid, cf, 'water', false);
            tmodel = TransportModel(model);
            modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
        else
            model = OverallCompositionCompositionalModel(G, rock, fluid, cf, 'water', false);
            modelFV = getSequentialModelFromFI(model);
        end
    else
        if opt.useGeneric
            model  = GenericNaturalVariablesModel(G, rock, fluid, cf, 'water', false);
            pmodel = PressureNaturalVariablesModel(G, rock, fluid, cf, 'water', false);
            tmodel = TransportModel(model);
            modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
        else
            model = NaturalVariablesCompositionalModel(G, rock, fluid, cf, 'water', false);
            modelFV = getSequentialModelFromFI(model);
        end
    end
    
    modelDG = cell(numel(opt.degree), 1);
    for dNo = 1:numel(opt.degree)
        tmodelDG = TransportModelDG(model, 'formulation', 'totalSaturation' , ...
                                           'degree'     , opt.degree(dNo), ...
                                           discArgs{:});
%         if 1
%             tmodelDG = TransportModelCompositionalDG(model, 'formulation', 'totalSaturation' , ...
%                                                'degree'     , opt.degree(dNo), ...
%                                                discArgs{:});
%         end
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
    
    setup = packSetup(state0, schedule, {{model}}, {{modelFV}}, {modelDG});
    
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function setup = spe1(args) %#ok

    opt = struct('degree', {{0, 1}}, 'k', {{[]}}, 'useGenericFV', true, 'ijk', [Inf, Inf, Inf]);
    [opt, discArgs] = merge_options(opt, args{:});
    
     if isempty(opt.k{1})
        k = cell(numel(opt.degree),1);
        [k{:}] = deal([]);
        opt.k = k;
    else
        opt.degree = cellfun(@(k) max(sum(k,2)), opt.k, 'UniformOutput', false);
    end
    
    [G, rock, fluid, deck, state0] = setupSPE1();
    gravity reset off
    
    G = computeGeometry(G);
    G = computeCellDimensions(G);
    G.cells.equal = false;
    
    model = selectModelFromDeck(G, rock, fluid, deck);
%     fluid = restrictRelperms(fluid);
    schedule = convertDeckScheduleToMRST(model, deck);
    if any(opt.ijk < Inf)
        pv0 = sum(poreVolume(G, rock));
        [ii, jj, kk] = gridLogicalIndices(G);
        opt.ijk = min([opt.ijk; G.cartDims], [],  1);
        keep = ii <= opt.ijk(1) & jj <= opt.ijk(2) & kk <= opt.ijk(3);
%         keep = ii <= 3 & jj <= 3 & kk == 1;
        G = extractSubgrid(G, keep);
        rock = extractSubrock(rock, keep);
        G = computeGeometry(G);
        G = computeCellDimensions(G);
        G.cells.equal = false;
        model = selectModelFromDeck(G, rock, fluid, deck);
    end
    pmodel = PressureModel(model);
    if opt.useGenericFV    
        tmodel = TransportModel(model);
        tmodel.formulation = 'missingPhase';
        tmodel.parentModel.useCNVConvergence = false;
        tmodel.parentModel.nonlinearTolerance = 1e-3;
        modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
    else
        m       = ThreePhaseBlackOilModel(G, rock, fluid, 'vapoil', model.vapoil, 'disgas', model.disgas);
        modelFV = getSequentialModelFromFI(m);
    end
        

    nls = NonLinearSolver('useLineSearch', true, 'enforceResidualDecrease', true, 'continueOnFailure', true, 'errorOnFailure', false);
    modelDG = cell(numel(opt.degree), 1);
    for dNo = 1:numel(opt.degree)                           
        tmodelDG     = TransportModelDG(model, ...
            'degree', opt.degree{dNo}, 'k', opt.k{dNo}, discArgs{:});
        tmodelDG.parentModel.drsMaxAbs = 200;
        tmodelDG.parentModel.useCNVConvergence = false;
        tmodelDG.parentModel.nonlinearTolerance = 1e-3;
        tmodelDG.parentModel.OutputStateFunctions = {};
        tmodelDG.formulation = 'missingPhase';
        modelDG{dNo} = SequentialPressureTransportModel(pmodel, tmodelDG, 'parentModel', model);
%         modelDG{dNo}.transportNonLinearSolver = nls;
    end
    
    % Convert the deck schedule into a MRST schedule by parsing the wells
    if any(opt.ijk < Inf)
        pv = sum(poreVolume(G, rock));
        schedule.control(1).W(2).cells = G.cells.num;
        schedule.control(1).W(1).val = schedule.control(1).W(1).val.*pv/pv0;
        schedule.control(1).W(2).val = schedule.control(1).W(2).val.*pv/pv0;
        state0.pressure = state0.pressure(G.cells.indexMap);
        state0.s        = state0.s(G.cells.indexMap,:);
        state0.rs       = state0.rs(G.cells.indexMap);
    end
    
    if 0
%         ix = 1:11;
        dt = 5*day;
        schedule.step.val = rampupTimesteps(sum(schedule.step.val), dt);
        schedule.step.control = ones(numel(schedule.step.val),1);
%         schedule.step.val     = schedule.step.val(ix);
%         schedule.step.control = schedule.step.control(ix);
    end
    setup = packSetup(state0, schedule, {{model}}, {{modelFV}}, {modelDG});
    
end

%-------------------------------------------------------------------------%
function setup = spe9(args) %#ok

    opt = struct('degree', {{0, 1}}, 'k', {{[]}}, 'useGenericFV', true, 'ijk', [Inf, Inf, Inf]);
    [opt, discArgs] = merge_options(opt, args{:});
    
     if isempty(opt.k{1})
        k = cell(numel(opt.degree),1);
        [k{:}] = deal([]);
        opt.k = k;
    else
        opt.degree = cellfun(@(k) max(sum(k,2)), opt.k, 'UniformOutput', false);
    end
    
    [G, rock, fluid, deck, state0] = setupSPE9();
    gravity reset off
    
    G = computeGeometry(G);
    G = computeCellDimensions2(G);
    G.cells.equal = false;
    
    model = selectModelFromDeck(G, rock, fluid, deck);
%     fluid = restrictRelperms(fluid);
    schedule = convertDeckScheduleToMRST(model, deck);
    if any(opt.ijk < Inf)
        [ii, jj, kk] = gridLogicalIndices(G);
        opt.ijk = min([opt.ijk; G.cartDims], [],  1);
        keep = ii >= G.cartDims(1) - opt.ijk(1) & jj >= G.cartDims(2) - opt.ijk(2) & kk <= opt.ijk(3);
%         keep = ii <= 3 & jj <= 3 & kk == 1;
        map = nan(G.cells.num,1);
        G = extractSubgrid(G, keep);
        map(G.cells.indexMap) = 1:G.cells.num;
        rock = extractSubrock(rock, keep);
        G = computeGeometry(G);
        G = computeCellDimensions2(G);
        G.cells.equal = false;
        model = selectModelFromDeck(G, rock, fluid, deck);
    end
    pmodel = PressureModel(model);
    if opt.useGenericFV    
        tmodel = TransportModel(model);
        tmodel.formulation = 'missingPhase';
        tmodel.parentModel.useCNVConvergence = false;
        tmodel.parentModel.nonlinearTolerance = 1e-3;
        modelFV = SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model);
    else
        m       = ThreePhaseBlackOilModel(G, rock, fluid, 'vapoil', model.vapoil, 'disgas', model.disgas);
        modelFV = getSequentialModelFromFI(m);
    end
        

    nls = NonLinearSolver('useLineSearch', true, 'enforceResidualDecrease', true, 'continueOnFailure', true, 'errorOnFailure', false);
    modelDG = cell(numel(opt.degree), 1);
    for dNo = 1:numel(opt.degree)
        disc         = DGDiscretization(modelFV, ...
                                   'degree', opt.degree{dNo}, 'k', opt.k{dNo}, discArgs{:});
        tmodelDG     = TransportModelDG(model, 'disc', disc);
        tmodelDG.parentModel.drsMaxAbs = 200;
        tmodelDG.parentModel.useCNVConvergence = false;
        tmodelDG.parentModel.nonlinearTolerance = 1e-3;
        tmodelDG.parentModel.OutputStateFunctions = {};
        tmodelDG.formulation = 'missingPhase';
        modelDG{dNo} = SequentialPressureTransportModel(pmodel, tmodelDG, 'parentModel', model);
%         modelDG{dNo}.transportNonLinearSolver = nls;
    end
    
    % Convert the deck schedule into a MRST schedule by parsing the wells
    if any(opt.ijk < Inf)
        pv = sum(poreVolume(G, rock));
        ctrl = schedule.control;
        for cNo = 1:numel(ctrl)
            isp = vertcat(ctrl(cNo).W.sign) < 0;
            prate0 = sum(vertcat(ctrl(cNo).W(isp).val));
            W = ctrl(cNo).W;
            keep = false(numel(W),1);
            for wNo = 1:numel(W)
                cells = map(W(wNo).cells);
                cells(isnan(cells)) = [];
                if ~isempty(cells)
                    keep(wNo) = true;
                    W(wNo).cells = cells;
%                     if any(strcmpi(W(wNo).type, {'rate', 'orat'}))
%                         W(wNo).val = W(wNo).val.*pv/pv0;
%                     end
                end
            end
            W = W(keep);
            ctrl(cNo).W = W;
            isp = vertcat(ctrl(cNo).W.sign) < 0;
            prate = sum(vertcat(ctrl(cNo).W(isp).val));
            isi = vertcat(ctrl(cNo).W.sign) > 0;
            ctrl(cNo).W(isi).val = ctrl(cNo).W(isi).val.*prate/prate0;
        end
        schedule.control = ctrl;
        state0.pressure = state0.pressure(G.cells.indexMap);
        state0.s        = state0.s(G.cells.indexMap,:);
        state0.rs       = state0.rs(G.cells.indexMap);
    end
    
    if 0
%         ix = 1:11;
        dt = 5*day;
        schedule.step.val = rampupTimesteps(sum(schedule.step.val), dt);
        schedule.step.control = ones(numel(schedule.step.val),1);
%         schedule.step.val     = schedule.step.val(ix);
%         schedule.step.control = schedule.step.control(ix);
    end
    setup = packSetup(state0, schedule, {{model}}, {{modelFV}}, {modelDG});
    
end

%-------------------------------------------------------------------------%
function setup = packSetup(state0, schedule, modelFI, modelFV, modelDG, varargin)

    setup = struct('state0'      , state0  , ...
                   'schedule'    , schedule, ...
                   'modelFI'     , modelFI , ...
                   'modelFV'     , modelFV , ...
                   'modelDG'     , modelDG , ...
                   'plot1d'      , false   );
    setup = merge_options(setup, varargin{:});
    
end

%-------------------------------------------------------------------------%
function fluid = restrictRelperms(fluid)
    names   = fieldnames(fluid)';
    for n = names
        if strcmpi(n{1}(1:2), 'kr')
            fluid.(n{1}) = @(s) fluid.(n{1})(min(max(s,0),1));
        end
    end
end
