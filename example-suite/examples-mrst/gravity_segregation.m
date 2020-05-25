function  [description, state0, model, schedule, options, plotOptions] = gravity_segregation(varargin)
    description = 'Gravity segregation example';
    options = struct('phases'   , 'WOG'       , ...
                     'ncells'   , 21          , ...
                     'nstep'    , 100         , ...
                     'type'     , 'horizontal', ...
                     'useRampup', true        );
    options = merge_options(options, varargin{:});
    description = [description, ' with ', options.phases                  , ...
                   ' and initially ', options.type, ' fluid phase contact'];  
    if nargout == 1, return; end
    nph = numel(options.phases);
    options.ncells = options.ncells - rem(options.ncells, nph);
    L = 100*meter;
    switch options.type
        case 'vertical'
            cartDims   = [options.ncells, 1, options.ncells];
            physDims   = [L, 1, L];
            pbar       = [1,1,1];
            sz         = [500, 500];
            view       = [0,0];
        case 'horizontal'
            cartDims   = [1, 1, options.ncells];
            physDims   = [1, 1, L];
            pbar       = [0.2,0.2,1];
            sz         = [300, 500];
            view       = [90,0];
    end
    
    gravity reset on
    G    = cartGrid(cartDims, physDims);
    G    = computeGeometry(G);
    rock = makeRock(G, 1*darcy, 0.3);

    phases    = 'WOG';
    fluidArgs = {'phases', phases                            , ...
                 'mu'    , [1,5,1]*centi*poise               , ...
                 'rho'   , [1500, 1000, 500]*kilogram/meter^3, ...
                 'n'     , [2,2,2]                           };
    active = ismember(phases, options.phases);
    fluidArgs(2:2:end) = cellfun(@(arg) arg(active), fluidArgs(2:2:end), 'UniformOutput', false);
    fluid = initSimpleADIFluid(fluidArgs{:}, 'cr', 1e-8);

    model = GenericBlackOilModel(G, rock, fluid, 'water', active(1), ...
                                                 'oil'  , active(2), ...
                                                 'gas'  , active(3));

    d = 1/nph;
    
    time = 5*year;
    dt   = rampupTimesteps(time, time/options.nstep, 5*options.useRampup);

    schedule = simpleSchedule(dt);

    [ii, ~, kk] = gridLogicalIndices(G);
    switch options.type
        case 'vertical'
            ll = ii;
        case 'horizontal'
            ll = kk;
    end

    state0 = initResSol(G, 1*atm, zeros(1,nph));
    for i = 1:nph
        cells = ll > max(ll)*d*(i-1) & ll <= max(ll)*d*i;
        sat = zeros(1,nph);
        sat(i) = 1;
        state0.s(cells,:) = repmat(sat, nnz(cells), 1);
    end
    % Plotting
    plotOptions = {'PlotBoxAspectRatio', pbar         , ...
                  'Projection'        , 'orthographic', ...
                  'View'              , view          , ...
                  'Size'              , sz            };
end