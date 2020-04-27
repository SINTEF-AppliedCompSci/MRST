function partition = getTopologicalFluxPartition(varargin)
    require matlab_bgl
    opt = struct('coarseGrid'      , []  , ...
                 'blockSize'       , -inf, ...
                 'processPartition', false);
    opt = merge_options(opt, varargin{:});
    partition = @(varargin) get(varargin{:}, opt);
end

function partition = get(model, state, state0, dt, drivingForces, opt)
    rmodel = model;
    if isa(model, 'WrapperModel')
        rmodel = rmodel.getReservoirModel();
    end
    if nargin == 2
        if ~isempty(state.coarseGrid)
            partition = state.coarseGrid.partition;
        else
            partition = ones(rmodel.G.cells.num,1);
        end
        return
    end
    v = state.flux;
    if ~isempty(opt.coarseGrid)
        G = opt.coarseGrid;
        W = [];
        cfsign = fineToCoarseSign(G);
        cfacesno = rldecode(1:G.faces.num, diff(G.faces.connPos), 2) .';
        vc = zeros(G.faces.num, size(state.flux, 2));
        
        for i = 1:size(v, 2)
            vc(:, i)   = accumarray(cfacesno, v(G.faces.fconn, i) .* cfsign);
        end
        v = vc;
        v = v(all(G.faces.neighbors > 0,2),:);
    else
        G = model.G;
        W = drivingForces.W;
        v = v(rmodel.operators.internalConn,:);
    end
    order = getTopologicalPermutation(G, v, 'W', W);
    if opt.blockSize < 0
        opt.blockSize = model.G.cells.num/20;
    end
    partition = sum(order >= 1:opt.blockSize:G.cells.num,2);
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
    if ~isempty(opt.coarseGrid)
        partition = partition(G.partition);
    end
end