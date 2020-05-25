function [description, state0, model, schedule, options, plotOptions] = qfs_wo(varargin)
    description = 'Quarter five-spot example on Cartesian grid';
    if nargout == 1, return; end
    % Quarter five-spot example on Cartesian grid
    options = struct('n', 50, 'coarseDims', [2,2], 'nkr', 2);
    options = merge_options(options, varargin{:});
    % Model
    G     = computeGeometry(cartGrid([1,1]*options.n, [1000, 1000]*meter));
    rock  = makeRock(G, 100*milli*darcy, 0.4);
    fluid = initSimpleADIFluid('phases', 'WO', 'n', [1,1]*options.nkr, 'mu', [1,1]*centi*poise, 'rho', [1,1]);
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    % Wells
    time = 2*year;
    rate = sum(poreVolume(G, rock))/time;
    bhp  = 50*barsa;
    W = [];
    W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate, 'comp_i', [1,0]);
    W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', bhp, 'comp_i', [1,0]);
    % Schedule
    schedule = simpleSchedule(rampupTimesteps(time, 30*day), 'W', W);
    % Initial state
    state0      = initResSol(G, 10*barsa, [0,1]);
    % Default plotting
    plotOptions = {};
end