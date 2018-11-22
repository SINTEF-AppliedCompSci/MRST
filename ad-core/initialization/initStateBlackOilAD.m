function [state, pressures] = initStateBlackOilAD(model, regions, varargin)
    opt = struct('pressure', []);
    opt = merge_options(opt, varargin{:});
    
    if ~iscell(regions)
        regions = {regions};
    end
    
    [rs, rv] = deal(0);
    vapoil = isprop(model, 'vapoil') && model.vapoil;
    disgas = isprop(model, 'disgas') && model.disgas;

    G = model.G;
    if disgas
        rs = zeros(G.cells.num, 1);
    end
    if disgas
        rv = zeros(G.cells.num, 1);
    end
    nph = sum(model.getActivePhases());
    state = struct('pressure', zeros(G.cells.num, 1), 'rs', rs, 'rv', rv, 's', zeros(G.cells.num, nph));
    
    watIx = model.getPhaseIndex('W');
    oilIx = model.getPhaseIndex('O');
    gasIx = model.getPhaseIndex('G');

    pressures = zeros(G.cells.num, nph);
    touched = false(G.cells.num, 1);
    for regNo = 1:numel(regions)
        region = regions{regNo};
        if ischar(region.cells)
            assert(numel(regions) == 1)
            region.cells = (1:model.G.cells.num)';
        end
        cells = region.cells;
        assert(~any(touched(cells)), 'Multiple regions defined in same cells.');
        touched(cells) = true;
        
        if isempty(opt.pressure)
            p = initializeEquilibriumPressures(model, region);
        else
            p = opt.pressure(cells, :);
        end
        z = model.G.cells.centroids(cells, 3);
        
        s = initializeEquilibriumSaturations(model, region, p);
        state.s(cells, :) = s;
        pressures(cells, :) = p;
        
        % Evalaute rel. perm.
        sat = cell(1, nph);
        pc = cell(1, nph);
        for i = 1:nph
            sat{i} = state.s(cells, i);
            pc{i} = region.pc_sign(i)*region.pc_functions{i}(state.s(cells, i));
        end
        kr = cell(1, nph);
        [kr{:}] = model.evaluateRelPerm(sat, 'cellInx', cells);
        kr = [kr{:}];
        numberOfMobile = sum(kr > 0, 2);
        maxSat = max(s, [], 2);
        singlePhaseMobile = numberOfMobile <= 1;
        
        toOil = true(size(p, 1), 1);
        if model.gas
            % If only gas is mobile, set oil pressure to the gas hydrostatic 
            % pressure minus the capillary pressure
            onlyGas = (kr(:, gasIx) > 0 & singlePhaseMobile) |...
                       (s(:, gasIx) == maxSat & numberOfMobile == 0);

            toOil(onlyGas) = false;
            state.pressure(cells(onlyGas)) = p(onlyGas, gasIx) - pc{gasIx}(onlyGas);
            if disgas
                po = p(:, oilIx);
                rsMax = model.fluid.rsSat(po, 'cellInx', cells);
                rs = region.rs(po, z);
                rs(rs > rsMax) = rsMax(rs > rsMax);
                state.rs(cells) = rs;
            end
        end
        if model.oil
            if vapoil
                pg = p(:, gasIx);
                rvMax = model.fluid.rvSat(po, 'cellInx', cells);
                rv = region.rv(pg, z);
                rv(rv > rvMax) = rvMax(rv > rvMax);
                state.rv(cells) = rv;
            end
        end
        if model.water
            % onlyWat = state.s(cells, watIx) == 1;
            onlyWat = (kr(:, watIx) > 0 & singlePhaseMobile) | ...
                      (s(:, watIx) == maxSat & numberOfMobile == 0);
            toOil(onlyWat) = false;
            state.pressure(cells(onlyWat)) = p(onlyWat, watIx) - pc{watIx}(onlyWat);
        end
        if model.oil
            state.pressure(cells(toOil)) = p(toOil, oilIx);
        end
    end
    if ~all(touched)
        warning('Regions did not cover all cells. Model only partially initialized.');
    end
end
