function [state0, model, schedule, keep] = removeDisconnectedCellGroups(state0, model, schedule)
% Remove cells which are part of disconnected parts of a domain, and not
% connected to any well.

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

    G = model.G;
    N = model.operators.N;
    W = schedule.control(1).W;
    for i = 1:numel(W)
        w = W(i);
        c = w.cells;
        if numel(c) > 1
            pairs = nchoosek(c, 2);
            N = [N; pairs]; %#ok
        end
    end
    require matlab_bgl
    A = getConnectivityMatrix(N);
    c = components(A);
    counts = accumarray(c, 1);
    fprintf('%d groups. Min: %d Max: %d (average: %1.0f cells)\n', ...
                        max(c), min(counts), max(counts), mean(counts));
    wellcells = vertcat(W.cells);
    cw = unique(c(wellcells));
    active_cells = ismember(c, cw);
    n = sum(active_cells);
    fprintf('%d groups are connected to wells (%d cells, %2.1f%%)\n', numel(cw), n, 100*n/G.cells.num);
    
    keep = false(G.cells.num, 1);
    keep(active_cells) = true;
    nc = sum(keep);
    
    renum = zeros(G.cells.num, 1);
    renum(active_cells) = (1:nc)';
    
    state0 = renum_state0(state0, keep, renum);
    model = renum_model(model, keep, renum);
    schedule = renum_schedule(schedule, keep, renum);
end

function state0 = renum_state0(state0, keep, renum)
    state0 = cellreduce(state0, keep);
    if isfield(state0, 'flux')
        state0 = rmfield(state0, 'flux');
    end
end

function model = renum_model(model, keep, renum)
    nc = sum(keep);
    N = model.operators.N;
    active_conn = all(keep(N), 2);
    N = renum(N(active_conn, :));
    T = model.operators.T(active_conn);
    
    pv = model.operators.pv(keep);
    
    model.operators.T_all = N;
    model.operators.internalConn = true(size(N, 1), 1);
    
    model.operators.pv = pv;
    model.operators.N = N;
    model.G.cells.centroids = model.G.cells.centroids(keep, :);
    model.G.cells.num = nc;
    if isfield(model.G.cells, 'indexMap')
        model.G.cells.indexMap = model.G.cells.indexMap(keep);
    end
    
    model.rock = cellreduce(model.rock, keep);
    
    model.operators = setupOperatorsTPFA(model.G, model.rock, ...
        'deck', model.inputdata, 'neighbors', N, 'trans', T, 'porv', pv);
    model.FlowPropertyFunctions = [];
    model.PVTPropertyFunctions = [];
    model.FlowDiscretization = [];
end

function schedule = renum_schedule(schedule, keep, renum)
    for i = 1:numel(schedule.control)
        W = schedule.control(i).W;
        for j = 1:numel(W)
            W(j).cells = renum(W(j).cells);
        end
        schedule.control(i).W = W;
    end
end

function x = cellreduce(x, keep)
    nc = numel(keep);
    flds = fieldnames(x);
    for i = 1:numel(flds)
        f = flds{i};
        if size(x.(f), 1) == nc
            x.(f) = x.(f)(keep, :);
        elseif isstruct(x.(f))
            x.(f) = cellreduce(x.(f), keep);
        end
    end
end
