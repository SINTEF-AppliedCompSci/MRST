function [state, pressures] = initStateBlackOilAD(model, regions, varargin)
    opt = struct('pressure', []);
    opt = merge_options(opt, varargin{:});
    
    [rs, rv] = deal(0);
    G = model.G;
    if model.disgas
        rs = zeros(G.cells.num, 1);
    end
    if model.vapoil
        rv = zeros(G.cells.num, 1);
    end
    nph = sum(model.getActivePhases());
    state = struct('pressure', zeros(G.cells.num, 1), 'rs', rs, 'rv', rv, 's', zeros(G.cells.num, nph));
    
    oilIx = model.getPhaseIndex('O');
    gasIx = model.getPhaseIndex('G');
    
    pressures = zeros(G.cells.num, nph);
    touched = false(G.cells.num, 1);
    for regNo = 1:numel(regions)
        region = regions{regNo};
        cells = region.cells;
        
        assert(~any(touched(cells)), 'Multiple regions defined in same cells.');
        touched(cells) = true;
        
        if isempty(opt.pressure)
            p = initializeEquilibriumPressures(model, region);
        else
            p = opt.pressure(cells, :);
        end
        
        if model.disgas
            rs = region.rs(p(:, oilIx), G.cells.centroids(cells, 3));
            state.rs(cells) = rs;
        end
        
        if model.vapoil
            rv = region.rv(p(:, gasIx), G.cells.centroids(cells, 3));
            state.rv(cells) = rv;
        end
        s = initializeEquilibriumSaturations(model, region, p);
        
        pressures(cells, :) = p;
        state.pressure(cells) = p(:, 2);
        state.s(cells, :) = s;
    end
    if ~all(touched)
        warning('Regions did not cover all cells. Model only partially initialized.');
    end
end