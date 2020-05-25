function [description, state0, model, schedule, options, plotOptions] = ifs_peaks_wo(varargin)
    description ...
        = ['Inverted five-spot pattern with perm/poro made from ', ...
           'repeated pattern of peaks'                             ];
    if nargout == 1, return; end
    options = struct('n', 51, 'tiles', 1, 'nkr', 2);
    options = merge_options(options, varargin{:});
    n = options.n - rem(options.n,2) + 1;
    % Model
    mrstModule add coarsegrid
    partition = partitionCartGrid([1,1]*n, [2,2]);
    m   = floor(n/2)+1;
    partition(m:n:end) = 5;
    partition(n*(m-1):n*m) = 5;
    k = zeros(n*n,1);
    for i = 1:4
        ix = find(partition == i);
        [xl, yl] = ind2sub([n,n], ix);
        XL = [xl, yl];
        XL = (XL - mean(XL))./(n/2)*6;
        if i == 2 || i == 4
            XL(:,1) = -XL(:,1);
        end
        if i > 2
            XL(:,2) = -XL(:,2);
        end
        k(ix) = peaks(XL(:,1), XL(:,2));
    end
    k = k./max(k);
    K = zeros(n+1);
    K(1:n,1:n) = reshape(k, n, n);
    K = repmat(K, options.tiles*[1,1]);
    K = reshape(K(1:end-1,1:end-1),[],1);
    N = n.*options.tiles+options.tiles-1;
    G = computeGeometry(cartGrid(N*[1,1], [1000, 1000]*meter));
    logperm = log10(100*milli*darcy) + K*3;
    poro    = 0.4*(1 + K*0.9);
    rock = makeRock(G, 10.^logperm, poro);
    fluid = initSimpleADIFluid('phases', 'WO'                         , ...
                               'n'     , [1,1].*options.nkr               , ...
                               'mu'    , [1,1]*centi*poise            , ...
                               'rho'   , [1000,1000]*kilogram/(meter^3));
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    % Wells
    time = 2*year;
    rate = sum(poreVolume(G, rock))/time;
    bhp  = 50*barsa;
    injector = @(W,i,j) verticalWell(W, G, rock, i, j, [], ...
                                         'type'  , 'rate', ...
                                         'val'   , rate  , ...
                                         'comp_i', [1,0] );
    producer = @(W,i,j) verticalWell(W, G, rock, i, j, [], ...
                                          'type'  , 'bhp', ...
                                          'val'   , bhp  , ...
                                          'comp_i', [1,0]);
    M   = floor(N/2)+1;
    W = [];
    W = injector(W, M, M);
    W = producer(W, 1, 1);
    W = producer(W, N, 1);
    W = producer(W, N, N);
    W = producer(W, 1, N);
    % Schedule
    schedule = simpleSchedule(rampupTimesteps(time, 30*day), 'W', W);
    % Initial state
    state0    = initResSol(G, bhp, [0,1]);
    plotOptions = {};
end