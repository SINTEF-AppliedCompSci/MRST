function schedule = upscaleSchedule(model, schedule, varargin)
    opt = struct('wellUpscaleMethod', 'sum');
    opt = merge_options(opt, varargin{:});

    for i = 1:numel(schedule.control)
        W = schedule.control(i).W;
        W_coarse = [];
        for j = 1:numel(W)
            w = handleWell(model, W(j), opt);
            W_coarse = [W_coarse; w]; %#ok
        end
        schedule.control(i).W = W_coarse;
    end
end

function Wc = handleWell(model, W, opt)
    % Handle a single well
    p = model.G.partition;
    Wc = W;
    
    s = W.cstatus;
    pc = p(W.cells);
    pc = pc(s);
    
    % cells
    [Wc.cells, firstInd, newMap] = uniqueStable(pc);
    nc = numel(Wc.cells);
    
    counts = accumarray(newMap, 1);
    
    % Take first direction uncritically (it shouldn't really matter when
    % the well index has been upscaled).
    dr = W.dir(s);
    Wc.dir = dr(firstInd);
    
    % Upscale well index using harmonic average...
    switch lower(opt.wellUpscaleMethod)
        case 'sum'
            fn = @(WI, map, counts) accumarray(map, WI);
        case 'harmonic'
            fn = @(WI, map, counts) 1./(accumarray(map, 1./WI(s))./counts);
        case 'mean'
            fn = @(WI, map, counts)accumarray(map, WI)./counts;
            otherwise
        error('Unknown')
    end
    Wc.WI = fn(W.WI(s), newMap, counts);

    % dZ
    z = model.G.cells.centroids(Wc.cells, 3);
    Wc.dZ = z - W.refDepth;
    
    % cstatus
    Wc.cstatus = true(nc, 1);
    
    % Extract topology.
    mp = [0; newMap];
    newtopo = mp(W.topo + 1);    
    % eliminate redundant connections due to cell collapsing in coarser
    % model
    newtopo = sort(newtopo, 2);
    newtopo = uniqueStable(newtopo, 'rows');
    newtopo = newtopo(newtopo(:,1) ~= newtopo(:, 2), :);
    
    Wc.topo = newtopo;
    Wc.parentIndices = firstInd;
end

function bc = handleBC(model, bc)
    if isempty(bc)
        return
    end
    error('NotImplemented');
    % Coarse face -> fine face map
    coarseFaceIx = rldecode(1:CG.faces.num, diff(CG.faces.connPos), 2) .';

    isPressure = strcmpi(bc.type, 'pressure');
    
    cf = unique(coarseFaceIx(bc.face));
end

function src = handleSRC(model, src)
    if isempty(src)
        return
    end
    error('NotImplemented');
end
