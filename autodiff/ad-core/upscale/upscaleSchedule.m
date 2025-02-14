function schedule = upscaleSchedule(model, schedule, varargin)
%Upscale a schedule to a coarser model
%
% SYNOPSIS:
%   schedule = upscaleSchedule(model, schedule)
%
% REQUIRED PARAMETERS:
%   model      - The coarse model the schedule is to be converted to.
%
%   schedule   - Schedule to be upscaled.
%
% OPTIONAL PARAMETERS:
%   'wellUpscaleMethod' - Upscaling method for well indices:
%                         - 'recompute' (default): Compute well indices in coarse grid.
%                         - 'sum': Sum of fine-scale WI.
%                         - 'harmonic': Harmonic average of WI.
%                         - 'mean': Arithmetic average of WI.
%                         - 'peaceman': Transmissibility-weighted WI (NEW).
%
%   'bcUpscaleMethod'   - Interpolation strategy for boundary conditions:
%                         - 'linear' (default), 'idw', 'mean', 'nearest'
%                         - 'flux-preserving' (NEW): Sum fluxes, interpolate pressure.
%
% RETURNS:
%   schedule   - Upscaled schedule for the coarse model.

    opt = struct('wellUpscaleMethod', 'peaceman', ...
                 'bcUpscaleMethod',   'flux');
    opt = merge_options(opt, varargin{:});

    for i = 1:numel(schedule.control)
        if isfield(schedule.control(i), 'W')
            W = schedule.control(i).W;
            W_coarse = [];
            for j = 1:numel(W)
                w = handleWell(model, W(j), opt);
                W_coarse = [W_coarse; w];
            end
            schedule.control(i).W = W_coarse;
        end

        if isfield(schedule.control(i), 'bc')
            schedule.control(i).bc = handleBC(model, schedule.control(i).bc, opt);
        end

%         if isfield(schedule.control(i), 'src')
%             schedule.control(i).src = handleSRC(model, schedule.control(i).src, opt);
%         end
    end
end

function Wc = handleWell(model, W, opt)
    % Handle well upscaling, including Peaceman correction
    p = model.G.partition;
    Wc = W;
    
    pc = p(W.cells);
    [Wc.cells, firstInd, newMap] = uniqueStable(pc);
    nc = numel(Wc.cells);
    counts = accumarray(newMap, 1);
    
    Wc.dir = W.dir(firstInd);
    if numel(W.r) > 1
        Wc.r = W.r(firstInd);
    end

    % Well index upscaling
    if strcmpi(opt.wellUpscaleMethod, 'recompute')
        Wc.WI = computeWellIndex(model.G, model.rock, Wc.r, Wc.cells, 'Dir', Wc.dir);
    else
        switch lower(opt.wellUpscaleMethod)
            case 'sum'
                fn = @(WI, map, counts) accumarray(map, WI);
            case 'harmonic'
                fn = @(WI, map, counts) 1./(accumarray(map, 1./WI)./counts);
            case 'mean'
                fn = @(WI, map, counts) accumarray(map, WI)./counts;
            case 'peaceman'
                T = computeTrans(model.G, model.rock);
                WI_weighted = W.WI .* T(Wc.cells);
                fn = @(WI, map, counts) accumarray(map, WI_weighted)./accumarray(map, T(Wc.cells));
            otherwise
                error(['Unknown well upscale mode: "', opt.wellUpscaleMethod, '"'])
        end
        Wc.WI = fn(W.WI, newMap, counts);
    end

    if model.G.griddim > 2
        z = model.G.cells.centroids(Wc.cells, 3);
    else
        z = zeros(numel(Wc.cells), 1);
    end
    Wc.dZ = z - W.refDepth;

    Wc.cstatus = true(nc, 1);
    Wc.fperf = newMap;
    Wc.parentIndices = firstInd;
end

function bc_coarse = handleBC(model, bc, opt)
    bc_coarse = [];
    if isempty(bc)
        return
    end
    CG = model.G;
    G = CG.parent;
    
    connCoarse = rldecode(1:CG.faces.num, diff(CG.faces.connPos), 2).';
    isFaceBC = false(G.faces.num, 1);
    isFaceBC(bc.face) = true;
    
    coarseFaceNo = zeros(G.faces.num, 1);
    coarseFaceNo(CG.faces.fconn) = connCoarse;
    
    coarseFacesBC = unique(connCoarse(isFaceBC(CG.faces.fconn)));
    for i = 1:numel(coarseFacesBC)
        cf = coarseFacesBC(i);
        isCurrentCoarse = coarseFaceNo == cf;
        act = isCurrentCoarse(bc.face);
        faces = bc.face(act);
        areas = G.faces.areas(faces);
        
        types = bc.type(act);
        type = types{1};
        values = bc.value(act);
        
        sat = [];
        comp = [];
        assert(all(strcmpi(type, types)), 'Mixed BC types on coarse face.');

        switch lower(type)
            case 'flux'
                if strcmpi(opt.bcUpscaleMethod, 'flux-preserving')
                    val = sum(values);
                else
                    val = sum(values);
                end
                if ~isempty(bc.sat)
                    sat = sum(bsxfun(@times, bc.sat(act, :), values), 1)./sum(values);
                end
                if isfield(bc, 'components')
                    comp = sum(bsxfun(@times, bc.components(act, :), values), 1)./sum(values);
                end
            case 'pressure'
                if ~isempty(bc.sat)
                    sat = sum(bsxfun(@times, bc.sat(act, :), areas), 1)./sum(areas);
                end
                if isfield(bc, 'components')
                    comp = sum(bsxfun(@times, bc.components(act, :), areas), 1)./sum(areas);
                end
                
                if strcmpi(opt.bcUpscaleMethod, 'flux-preserving')
                    val = sum(bsxfun(@times, values, areas))/sum(areas);
                else
                    xq = CG.faces.centroids(cf, :);
                    x = G.faces.centroids(faces, :);
                    if numel(faces) == 1
                        val = values;
                    else
                        switch lower(opt.bcUpscaleMethod)
                            case 'idw'
                                val = interpolateIDW(x, values, xq, 2);
                            case 'mean'
                                val = mean(values.*areas)./sum(areas);
                            case 'nearest'
                                val = nearestNeighbor(x, values, xq);
                            otherwise
                                val = interp1(x(:), values(:), xq, 'linear', 'extrap');
                        end
                        if isempty(val) || any(~isfinite(val))
                            val = nearestNeighbor(x, values, xq);
                        end
                    end
                end
            otherwise
                error(['Unknown BC type: "', type, '"']);
        end
        
        bc_coarse = addBC(bc_coarse, cf, type, val, 'sat', sat);
        if isfield(bc, 'components')
            bc_coarse.components = [bc_coarse.components; comp];
        end
    end
end
