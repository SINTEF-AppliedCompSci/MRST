function regions = getInitializationRegionsDeck(model, deck)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    n = size(deck.SOLUTION.EQUIL, 1);
    model = model.validateModel();

    regions = cell(n, 1);
    for regionIx = 1:n
        % Get all combinations of PVT/SATNUM regions for this block
        if isfield(deck.REGIONS, 'EQLNUM')
            cells = find(deck.REGIONS.EQLNUM(model.G.cells.indexMap) == regionIx);
        else
            cells = (1:model.G.cells.num)';
        end
        eql = deck.SOLUTION.EQUIL(regionIx, :);
        
        sat = getSATNUM(model, deck, cells);
        pvt = getPVTNUM(model, deck, cells);
        
        pairs = unique([sat, pvt], 'rows');
        nreg_comb = size(pairs, 1);
        
        sub_regions = cell(nreg_comb, 1);
        for j = 1:nreg_comb
            satnum = pairs(j, 1);
            pvtnum = pairs(j, 2);
            subcells = cells(sat == satnum & pvt == pvtnum);
            
            sub_regions{j} = getRegion(model, deck, eql, subcells, regionIx, satnum, pvtnum);
        end
        regions{regionIx} = sub_regions;
    end
    % Flatten local regions
    regions = vertcat(regions{:});
end

function reg = getSATNUM(model, deck, cells)
    if isfield(deck.REGIONS, 'SATNUM')
        reg = deck.REGIONS.SATNUM(model.G.cells.indexMap(cells));
    else
        reg = ones(numel(cells), 1);
    end
end


function reg = getPVTNUM(model, deck, cells)
    if isfield(deck.REGIONS, 'PVTNUM')
        reg = deck.REGIONS.PVTNUM(model.G.cells.indexMap(cells));
    else
        reg = ones(numel(cells), 1);
    end
end

function region = getRegion(model, deck, eql, cells, regionIx, satnum, pvtnum)
    rs_method = eql(7);
    rv_method = eql(8);
    contacts = eql([3, 5]);
    contacts_pc = eql([4, 6]);
    p_datum = eql(2);
    if model.water && model.gas && ~model.oil
        act = [true, false, false];
    else
        act = [model.water & model.oil, model.oil & model.gas];
    end
    % Handle specifics
    hasVapoil = isprop(model, 'vapoil') && model.vapoil;
    hasDisgas = isprop(model, 'disgas') && model.disgas;
    hasEOS = isprop(model, 'EOSModel');
    [rs, rv] = deal([]);
    if hasDisgas
        if rs_method <= 0
            rsSat = model.fluid.rsSat;
            if iscell(rsSat)
                rs =  @(p, z) 0*p + model.fluid.rsSat{pvtnum}(p_datum);
            else
                rs =  @(p, z) 0*p + model.fluid.rsSat(p_datum);
            end
        else
            if isfield(deck.SOLUTION, 'RSVD')
                rsvd = deck.SOLUTION.RSVD{regionIx};
            else
                assert(isfield(deck.SOLUTION, 'PBVD'));
                rsSat = model.fluid.rsSat;
                if iscell(rsSat)
                    rsSat = rsSat{pvtnum};
                end
                pbvd = deck.SOLUTION.PBVD{regionIx};
                rsvd = [pbvd(:, 1), rsSat(pbvd(:, 2))];
            end
            F = griddedInterpolant(rsvd(:, 1), rsvd(:, 2), 'linear', 'nearest');
            rs = @(p, z) F(z);
        end
    end
    if hasVapoil
        if rv_method <= 0
            % Oil pressure at gas-oil contact + capillary pressure there
            pg_goc = p_datum + eql(6);
            rvSat = model.fluid.rvSat;
            if iscell(rvSat)
                rv =  @(p, z) 0*p + model.fluid.rvSat{pvtnum}(pg_goc);
            else
                rv =  @(p, z) 0*p + model.fluid.rvSat(pg_goc);
            end
        else
            assert(isfield(deck.SOLUTION, 'RVVD'));
            rvvd = deck.SOLUTION.RVVD{regionIx};
            F = griddedInterpolant(rvvd(:, 1), rvvd(:, 2), 'linear', 'nearest');
            rv = @(p, z) F(z);
        end
    end

    if hasEOS
        z_method = eql(10);
        assert(z_method == 1);
        zvd = deck.PROPS.ZMFVD{regionIx};
        z_fn = @(p, z) interpolateDepthTable(zvd(:, 1), zvd(:, 2:end), z);
        if isfield(deck.PROPS, 'RTEMP')
            T = deck.PROPS.RTEMP;
            T_fn = @(p, z) repmat(T, size(p));
        else
            tvd = deck.PROPS.TEMPVD{regionIx};
            T_fn = @(p, z) interpolateDepthTable(tvd(:, 1), tvd(:, 2), z);
        end
        contacts(2) = min(model.G.cells.centroids(:, 3)) - 10;

        region = getInitializationRegionsCompositional(model, contacts(act),...
            'cells', cells, 'datum_pressure', p_datum, ...
            'datum_depth', eql(1), 'contacts_pc', contacts_pc(act), ...
            'x', z_fn, 'y', z_fn, 'T', T_fn);
        region.T = T_fn;
        region.z = z_fn;
    else
        region = getInitializationRegionsBlackOil(model, contacts(act),...
            'cells', cells, 'datum_pressure', p_datum, ...
            'datum_depth', eql(1), 'contacts_pc', contacts_pc(act), ...
            'rs', rs, 'rv', rv);
    end
    if hasDisgas
        region.rs = rs;
    end
    
    if hasVapoil
        region.rv = rv;
    end
end

function fq = interpolateDepthTable(x, f, xq)
    nd = size(f, 2);
    n = size(xq, 1);
    fq = zeros(n, nd);
    for i = 1:nd
        T = griddedInterpolant(x, f(:, i), 'linear', 'nearest');
        fq(:, i) = T(xq);
    end
end
