function [state, pressures] = initStateBlackOilAD(model, regions, varargin)
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
    % Make sure that model is ready for evaluation of properties
    model = model.validateModel();
    opt = struct('pressure', []);
    opt = merge_options(opt, varargin{:});
    
    if ~iscell(regions)
        regions = {regions};
    end
    vapoil = isprop(model, 'vapoil') && model.vapoil;
    disgas = isprop(model, 'disgas') && model.disgas;
    compositional = isprop(model, 'EOSModel');

    G = model.G;
    nph = sum(model.getActivePhases());
    state = struct('pressure', zeros(G.cells.num, 1), 's', zeros(G.cells.num, nph));
    if disgas
        state.rs = zeros(G.cells.num, 1);
    end
    if vapoil
        state.rv = zeros(G.cells.num, 1);
    end
    if compositional
        ncomp = model.EOSModel.CompositionalMixture.getNumberOfComponents();
        state.components = zeros(G.cells.num, ncomp);
        state.T = zeros(G.cells.num, 1);
    end
    watIx = model.getPhaseIndex('W');
    oilIx = model.getPhaseIndex('O');
    gasIx = model.getPhaseIndex('G');
    % Reference phase should be oil if present, otherwise it we assume it
    % to be gas, then just the first phase
    if model.oil
        refIx = oilIx;
    elseif model.gas
        refIx = gasIx;
    else
        refIx = 1;
    end

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
        pc = cell(1, nph);
        for i = 1:nph
            pc{i} = region.pc_sign(i)*region.pc_functions{i}(state.s(cells, i));
        end
        m = model;
        m.FlowPropertyFunctions = model.FlowPropertyFunctions.subset(cells);

        tmp_state = struct('s', s, 'pressure', p(:, 1));
        tmp_state = model.initStateFunctionContainers(tmp_state);
        try
            kr = value(m.getProp(tmp_state, 'RelativePermeability'));
        catch
            warning(['Unable to evaluate relative permeability during ', ...
                    'initialization. The phase pressures might be off where', ...
                    ' saturations do not correspond to mobile phases']);
            kr = s;
        end
        maxSat = max(s, [], 2);
        referenceImmobile = kr(:, refIx) < 1e-8;
        
        toReferencePhase = true(size(p, 1), 1);
        if model.gas
            % If only gas is mobile, set oil pressure to the gas hydrostatic 
            % pressure minus the capillary pressure
            gasMajority = s(:, gasIx) == maxSat;
            onlyGas = gasMajority & referenceImmobile;

            toReferencePhase(onlyGas) = false;
            state.pressure(cells(onlyGas)) = p(onlyGas, gasIx) - pc{gasIx}(onlyGas);
            if disgas
                po = p(:, oilIx);
                if iscell(model.fluid.rsSat)
                    rsSatF = model.fluid.rsSat{region.pvt_region};
                else
                    rsSatF = model.fluid.rsSat;
                end
                rsMax = rsSatF(po);
                rs = region.rs(po, z);
                rs(rs > rsMax) = rsMax(rs > rsMax);
                sg = s(:, gasIx);
                rs(sg > 0) = rsMax(sg > 0);
                state.rs(cells) = rs;             
            end
        end
        if model.oil
            if vapoil
                pg = p(:, gasIx);
                if iscell(model.fluid.rvSat)
                    rvSatF = model.fluid.rvSat{region.pvt_region};
                else
                    rvSatF = model.fluid.rvSat;
                end
                rvMax = rvSatF(po);
                rv = region.rv(pg, z);
                rv(rv > rvMax) = rvMax(rv > rvMax);
                so = s(:, oilIx);
                rv(so > 0) = rvMax(so > 0);
                state.rv(cells) = rv;
            end
        end
        if model.water
            watMajority = s(:, watIx) == maxSat;
            onlyWat = watMajority & referenceImmobile;
            toReferencePhase(onlyWat) = false;
            state.pressure(cells(onlyWat)) = p(onlyWat, watIx) - pc{watIx}(onlyWat);
        end
        % Set remaining pressure to reference phase
        state.pressure(cells(toReferencePhase)) = p(toReferencePhase, refIx);
        
        if compositional
            if isfield(region, 'z')
                state.components(cells, :) = region.z(state.pressure(cells), z);
            end
            if isfield(region, 'T')
                state.T(cells, :) = region.T(state.pressure(cells), z);
            end
        end
    end
    if ~all(touched)
        warning('Regions did not cover all cells. Model only partially initialized.');
    end
    if isprop(model, 'polymer')
        state.cp = zeros(G.cells.num, 1);
    end
    if isprop(model, 'surfactant')
        state.cs = zeros(G.cells.num, 1);
    end
end
