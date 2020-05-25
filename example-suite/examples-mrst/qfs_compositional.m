function [description, state0, model, schedule, options, plotOptions] = qfs_compositional(varargin)
    description = ...
        ['Quarter five-spot with CO2 injection into methane, n-decane' ...
         'and CO2. Slightly modified from  from Klemetsdal et al,'     ...
         'SPE RSC 2019, doi: 10.2118/193934-ms'                        ];
    if nargout == 1, return; end
    
    options = struct('n', 60, 'nstep', 100);
    options = merge_options(options, varargin{:});
    mrstModule add reordering compositional
    gravity reset off
    % Model
    nx = options.n; ny = options.n; nz = 1;
    G = computeGeometry(cartGrid([nx, ny, nz], [1000, 1000, 10]));
    K = 0.1*darcy;
    rock = makeRock(G, K, 0.25);
    fluid = initSimpleADIFluid('n', [2, 3],...
                               'phases', 'OG', ...
                               'rho', [100, 100]);
    [cf, info] = getCompositionalFluidCase('simple');
    model = GenericOverallCompositionModel(G, rock, fluid, cf, 'water', false);
    % Schedule
    time = 5*year;
    irate = 100*sum(model.operators.pv)/time;
    W = [];
    W = addWell(W, G, rock, 1, 'type', 'rate', 'val', irate, 'comp_i', [1, 0]);
    W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 25*barsa, 'comp_i', [1, 0]);
    for i = 1:numel(W)
        W(i).components = info.injection;
    end
    dt = rampupTimesteps(time, time/options.nstep);
    schedule = simpleSchedule(dt, 'W', W);
    state0 = initCompositionalState(G, 50*barsa, info.temp, [1, 0], info.initial, model.EOSModel);
    % Plotting
    plotOptions = {'PlotBoxAspectRatio', [1,1,1], 'View', [0, 90], 'Projection', 'orthographic'};
end