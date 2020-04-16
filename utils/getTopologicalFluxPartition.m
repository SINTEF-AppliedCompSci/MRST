function partition = getTopologicalFluxPartition(varargin)
    require matlab_bgl
    opt = struct('blockSize'       , -inf, ...
                 'processPartition', false);
    opt = merge_options(opt, varargin{:});
    if opt.blockSize < 0
        opt.blockSize = model.G.cells.num/20;
    end
    partition = @(varargin) get(varargin{:}, opt);
end

function partition = get(model, state, state0, dt, drivingForces, opt)
    rmodel = model;
    if isa(model, 'WrapperModel')
        rmodel = rmodel.getReservoirModel();
    end
    if nargin == 2
        partition = ones(rmodel.G.cells.num,1);
        return
    end
    v = state.flux(rmodel.operators.internalConn,:);
    order = getTopologicalPermutation(model.G, v, 'W', drivingForces.W);
    partition = sum(order >= 1:opt.blockSize:model.G.cells.num,2);
    W = drivingForces.W;
    if ~isempty(W)
        M = getConnectivityMatrix(rmodel.operators.N);
        for i = 1:numel(W)
            cells = false(rmodel.G.cells.num,1);
            cells(W(i).cells) = true;
            cell0 = find(cells);
            cells = cells | M*cells;
            partition(cells) = partition(cell0(1));
        end
    end
    partition = compressPartition(partition);
    if opt.processPartition
        partition = processPartition(rmodel.G, partition);
    end
end