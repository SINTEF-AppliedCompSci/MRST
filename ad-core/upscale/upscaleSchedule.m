function schedule = upscaleSchedule(model, schedule, varargin)
    

    for i = 1:numel(schedule.control)
        W = schedule.control(i).W;
        for j = 1:numel(W)
            W(j) = handleWell(model, W(j));
        end
        schedule.control(i).W = W;
    end
end

function Wc = handleWell(model, W)
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
    Wc.WI = 1./(accumarray(newMap, 1./W.WI(s))./counts);

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
    newtopo = uniqueStable(newtopo, 'rows');
    newtopo = newtopo(newtopo(:,1) ~= newtopo(:, 2), :);
    
    Wc.topo = newtopo;
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
