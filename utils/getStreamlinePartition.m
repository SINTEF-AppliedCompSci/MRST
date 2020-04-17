function partition = getStreamlinePartition(varargin)
    
    require matlab_bgl
    opt = struct('state', [], 'W', []);
    opt = merge_options(opt, varargin{:});
    partition = @(varargin) get(varargin{:}, opt);

end

function partition = get(model, state, state0, dt, drivingForces, opt)
    rmodel = model;
    if isa(model, 'WrapperModel')
        rmodel = rmodel.getReservoirModel();
    end
    if nargin == 2
        opt = state;
        state = opt.state;
    end
    G = rmodel.G;
    cells = (G.cartDims(1):G.cartDims(1)-1:prod(G.cartDims)-1)';
    cells = cells(5:10:end);
    
    SF = pollock(G, state, cells);
    SR = pollock(G, state, cells, 'reverse', true);
    
    S = cellfun(@(sf, sr) vertcat(sf, sr), SF, SR, 'UniformOutput', false);
    
    x = G.cells.centroids;
    d = nan(G.cells.num, numel(S));
    for i = 1:numel(S)
        d(:,i) = min(pdist2(x, S{i}),[],2);
    end
    [~, partition] = min(d, [], 2);
    if nargin == 2
        W = opt.W;
    else
        W = drivingForces.W;
    end
    if ~isempty(W)
        x = G.cells.centroids;
        rWell = 150;
%         M = getConnectivityMatrix(rmodel.operators.N);
        for i = 1:numel(W)
            cells = pdist2(x(W(i).cells,:), x) < rWell;
%             cells = false(rmodel.G.cells.num,1);
%             cells(W(i).cells) = true;
%             cell0 = find(cells);
%             for i = 1:10
%                 cells = cells | M*cells;
%             end
            partition(cells) = max(partition) + 1;
        end
    end
    partition = compressPartition(partition);
    partition = processPartition(rmodel.G, partition);
%     T = rmodel.operators.T_all;
%     partition = mergeBlocksByConnections(G, partition, T, 20);
end