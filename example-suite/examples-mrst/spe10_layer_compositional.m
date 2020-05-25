function [description, state0, model, schedule, options, plotOptions] = spe10_layer_compositional(varargin)    
    description ...
        = ['Water-alternating gas injection in layer 1 of SPE10 Model 2 ', ...
           'using a six-component model. Example from Moncorgé et al, '  , ...
           'J. Comput. Phys, 2018, doi: 10.1016/j.jcp.2018.05.048'       ];
    if nargout == 1, return; end
   
    options = struct('diagonalAD', false, 'dt', 25*day);
    options = merge_options(options, varargin{:});

    mrstModule add spe10 compositional deckformat ad-core ad-props

    fldr = fullfile(mrstDataDirectory(), 'arthur_2d');
%     fldr = fullfile(mrstPath('ddc'),'experimental', 'data', 'arthur_2d');
    fn = fullfile(fldr, 'TEST1.DATA');

    deck = readEclipseDeck(fn);
    deck = convertDeckUnits(deck);

    G = SPE10_setup(1);
    G = computeGeometry(G);

    rock = initEclipseRock(deck);
    rock = compressRock(rock, G.cells.indexMap);
    eos  = initDeckEOSModel(deck);

    fluid = initDeckADIFluid(deck);

    inj_comp =  [1.0 0.0 0.0 0.0 0.0 0.0];

    p_std = 14.7*psia;
    T_std = 288.7;
    mc = sum(inj_comp.*eos.fluid.molarMass);
    R = 8.3144598;
    rhoL = p_std/(R*T_std);
    rhoL = rhoL.*mc;
    rhoV = rhoL;

    fluid.rhoOS = rhoL;
    fluid.rhoGS = rhoV;

    model = GenericOverallCompositionModel(G, rock, fluid, eos.fluid, 'water', true);
    if options.diagonalAD
        model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
    end
    schedule = convertDeckScheduleToMRST(model, deck);

    nc = numel(schedule.control);
    schedule0 = schedule;

    dt = cell(nc, 1);
    for i = 1:nc
        t_loc = sum(schedule.step.val(schedule.step.control == i));
        dt{i} = rampupTimesteps(t_loc, options.dt, 8);
    end

    schedule.step.val = vertcat(dt{:});
    schedule.step.control = rldecode((1:nc)', cellfun(@numel, dt));
    schedule.s0 = schedule0;


    for i = 1:3
        schedule.control(i).W(1).components = inj_comp;
        schedule.control(i).W(2).components = inj_comp;
    end

    z0 = [0.5 0.03 0.07 0.2 0.15 0.05];
    T = 344;
    state0 = initCompositionalState(G, 4000*psia, T, [0, 0, 1], z0, eos);

    % Plotting
    plotOptions = {'PlotBoxAspectRatio', [1,1.83,0.3]  , ...
                   'Projection'        , 'orthographic', ...
                   'View'              , [0, 90]       , ...
                   'Size'              , [500, 800]    };
    
end