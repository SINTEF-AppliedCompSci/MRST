function [W, Gsub] = processWellboreTrajectories(G, rock, coords, varargin)
    
    % Process input arguments
    %---------------------------------------------------------------------%
    opt = struct('names'    , {{}}, ...
                 'radius'   , []  , ...
                 'refDepth' , []  , ...
                 'roughness', 1e-4);
    [opt, wellsGiven, varargin] = processInput(coords, opt, varargin{:});
    %---------------------------------------------------------------------%
    
    % Process trajectories
    %---------------------------------------------------------------------%
    require wellpaths
    if ~isfield(G.faces, 'bbox'), G = addBoundingBoxFields(G); end
    
    nw = numel(coords);
    W  = [];
    if wellsGiven, W = coords; coords = cell(nw, 1); end
    
    trajectories = cell(nw,1);
    for i = 1:nw
        if wellsGiven
            W(i) = topoSortWellCells(G, W(i)); %#ok
            coords{i} = G.cells.centroids(W(i).cells,:);
            traj = computeTraversedCells(G, coords{i});
            if isfield(G.cells, 'hybrid') && any(G.cells.hybrid(W(i).cells))
                traj = addHybridCells(G, W(i), traj);
            end
            trajectories{i} = traj;
        else
            args = {};
            if ~isempty(opt.names)   , args = [args, {'name'    , opt.names{i}}   ]; end %#ok
            if ~isempty(opt.radius)  , args = [args, {'radius'  , opt.radius(i)}  ]; end %#ok
            if ~isempty(opt.refDepth), args = [args, {'refDepth', opt.refDepth(i)}]; end %#ok
            [W, trajectories{i}] = addWellFromTrajectory(W, G, rock, coords{i}, args{:}, varargin{:});
        end
        trajectories{i}.roughness = opt.roughness(i);
    end
    [W.trajectory] = deal(trajectories{:});
    %---------------------------------------------------------------------%
    
    % Construct wellbore grid
    %---------------------------------------------------------------------%
    % Get traversed cell indices
    nc    = arrayfun(@(W) numel(W.cells), W);
    cells = vertcat(W.cells);
    Gsub  = extractSubgridLocal(G, cells); clear G;
    [~, order] = sort(cells);
    % Tag cells for each well
    wellNo = rldecode(1:numel(trajectories), nc, 2)';
    Gsub.cells.wellNo = wellNo(order);
    Gsub.cells.order = order;
    % Find top cell and top face for each well
    [Gsub.cells.topCell, Gsub.faces.topFace] = deal(nan(numel(trajectories),1));
    for i = 1:nw
        % Get top cell (well cells should be sorted at this point)
        topCell = find(Gsub.cells.global ==  W(i).cells(1));
        faces = Gsub.cells.faces(Gsub.cells.facePos(topCell):Gsub.cells.facePos(topCell+1)-1);
        if wellsGiven
            dim = 'xyz' == W(i).dir(1);
            cells = find(Gsub.cells.wellNo == i);
            dir = diff(Gsub.cells.centroids(cells(1:2),dim));
            if isempty(dir) || dir > 0
                fun = @min;
            else
                fun = @max;
            end
            [~, fix] = fun(Gsub.faces.centroids(faces,dim));
        else
            [~, fix] = min(minPdist2(Gsub.faces.centroids(faces,:), coords{i}(1,:)));
        end
        topFace = faces(fix(1));
        Gsub.faces.topFace(i) = topFace;
        Gsub.cells.topCell(i) = topCell;
    end
    %---------------------------------------------------------------------%

end

%-------------------------------------------------------------------------%
function [opt, wellsGiven, varargin] = processInput(coords, opt, varargin)

    nw = numel(coords);
    wellsGiven = isstruct(coords);
    
    if ~wellsGiven
        opt.radius   = 10*centi*meter;
        opt.refDepth = 0*meter;
        opt.names    = cellfun(@(n) sprintf('W%d', n), ...
                                num2cell(1:nw), 'UniformOutput', false);
    end
    [opt, varargin] = merge_options(opt, varargin{:});

    if numel(opt.radius)    == 1, opt.radius    = repmat(opt.radius   , nw, 1); end
    if numel(opt.refDepth)  == 1, opt.refDepth  = repmat(opt.refDepth , nw, 1); end
    if numel(opt.roughness) == 1, opt.roughness = repmat(opt.roughness, nw, 1); end
    assert(any(numel(opt.radius)    == [0,nw]), 'Please provide one radius per well'          );
    assert(any(numel(opt.refDepth)  == [0,nw]), 'Please provide one reference depth per well' );
    assert(any(numel(opt.roughness) == [0,nw]), 'Please provide one roughness factor per well');
    
end

%-------------------------------------------------------------------------%
function Gsub = extractSubgridLocal(G, cells)
    
    Gsub = extractSubgrid(G, cells);
    
    if isfield(G, 'hybridNeighbors')
        
        % Set hybrid filed
        Gsub.cells.hybrid = G.cells.hybrid(Gsub.cells.global);
        Gsub.faces.hybrid = G.faces.hybrid(Gsub.faces.global);
        
        Gsub.hybridNeighbors = G.hybridNeighbors;
        Gsub.hybridNeighbors.neighbors = [];
        Gsub.hybridNeighbors.facePos   = [];
        Gsub.hybridNeighbors.faces     = [];
        Gsub.hybridNeighbors.n         = [];
        
    end
    
end

%-------------------------------------------------------------------------%
function W = topoSortWellCells(G, W)
% Topologically sort well cells. Mostly for visualization - should not have
% any effect on drawdown computations.

    cells = W.cells;
    n = numel(cells);
    if n == 1, return; end
    
    cells0 = cells; %#ok For debugging
    
    loc2glob = cells;
    glob2loc = nan(G.cells.num, 1); glob2loc(cells) = (1:n)';

    facePos = mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1);
    faces   = G.cells.faces(facePos,1);
    N       = G.faces.neighbors(faces,:);
    keep = all(ismember(N, cells), 2);
    N = N(keep,:); 
    N = glob2loc(N);
    N = unique(N, 'rows');
    N = [N; fliplr(N)];
    A = sparse(N(:,1), N(:,2), 1, n, n);
   
    endpoints = find(sum(A,2) == 1);
    
    [~, ix] = min(G.cells.centroids(endpoints,3));
    first = endpoints(ix);
 
    cells = nan(n,1);
    current = false(n,1);
    current(first) = true;
    
    processed = false(n,1);
    pno = 1;
    while any(~processed)
        cells(pno) = find(current);
        processed(current) = true;
        candidates = A*current > 0 & ~processed;
        if nnz(candidates) > 1
            warning('Branching wells not supported')
            cix = find(candidates);
            candidates(cix(2:end)) = false;
            processed(cix(2:end)) = true;
        end
        
        current = candidates;
        pno = pno + 1;
    end
    cells = cells(~isnan(cells));
    order = cells;
    cells = loc2glob(cells);
    W.cells = cells;
    
    names = fieldnames(W);
    for name = names'
        v = W.(name{1});
        if size(v,1) == n
            W.(name{1}) = v(order,:);
        end
    end

end

function traj = addHybridCells(G, W, traj)

    hybrid = ~ismember(W.cells, traj.cell);
    nc = numel(W.cells);
    
    assert(nnz(~hybrid) == nnz(~G.cells.hybrid(W.cells)))
    traj.cell = W.cells;
    traj.vec = interp1(find(~hybrid), traj.vec, (1:nc)');
    traj.vec(hybrid,3) = G.cells.aperture(W.cells(hybrid));
    w0 = traj.weight;
    traj.weight = zeros(nc,1);
    traj.weight(~hybrid) = w0;
    traj.weight(hybrid) = 1;
    
end