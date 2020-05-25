function [description, state0, model, schedule, options, plotOptions] = simple_1d(varargin)
    description = 'Simple 1D Buckley Leverett displacement';
    if nargout == 1, return; end
    options = struct('n', 100, 'nkr', 2, 'cfl', 1);
    options = merge_options(options, varargin{:});
    description = [description, ' with CFL = ', num2str(options.cfl)];
    % Model
    G     = computeGeometry(cartGrid([options.n,1]));
    rock  = makeRock(G, 1, 1);
    fluid = initSimpleADIFluid('phases', 'WO', 'n', [1,1]*options.nkr, 'mu', [1,1], 'rho', [1,1]);
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    % Set up wells
    W = [];
    time = options.n/options.cfl;
    rate = sum(poreVolume(G, rock))/time;
    W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate, 'WI', 9999, 'comp_i', [1,0]);
    W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 0.5, 'WI', 9999, 'comp_i', [1,0]);
    % Make schedule
    schedule = simpleSchedule(rampupTimesteps(options.n/options.cfl, 1, 0), 'W', W);
    % Initial state
    state0 = initResSol(G, 1, [0,1]);
    % Plotting
    plotOptions = {'plot1d'            , true      , ...
                   'lockCaxis'         , true      , ...
                   'Size'              , [800, 350], ...
                   'PlotBoxAspectRatio', [2.5,1,1] , ...
                   'YLim', [0,1]                   };
end