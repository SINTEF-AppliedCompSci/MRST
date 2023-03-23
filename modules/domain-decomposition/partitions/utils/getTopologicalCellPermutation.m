function [order, M, cyclesPresent, cycleSize] = getTopologicalCellPermutation(G, v, varargin)
    % Get cell permutation for a grid G based on the topological order
    % induced by the edges v. Input v has one entry per internal
    % connection, and may have multiple columns. The sign of an edge is
    % interpreted in the same way as intercell fluxes.

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    require matlab_bgl
    opt = struct('W'             , []   , ... % Well structure
                 'padWells'      , true , ... % Pad wells with one cell
                 'cycleTolerance', -inf);     % Threshold for edge weight for cycle across single interface
                 
    opt = merge_options(opt, varargin{:});

    N = G.faces.neighbors;
    internalConn = all(N ~= 0, 2);
    N = N(internalConn, :);
    if isstruct(v)
        v = v.flux(internalConn,:);
    end
    nonzero = abs(v) > opt.cycleTolerance; % Disregard cycles if edge weights are very small
    pos = any(v > 0 & nonzero, 2);
    neg = any(v < 0 & nonzero, 2);
    countercurrent = pos & neg;
    
    I0 = [N(pos, 1); N(neg, 2)];
    J0 = [N(pos, 2); N(neg, 1)];
    
    if ~isempty(opt.W)
        % Wells were sent in
        C = getConnectivityMatrix(N);
        for i = 1:numel(opt.W)
            if opt.padWells
                cl = false(G.cells.num,1);
                cl(opt.W(i).cells) = true;
                cl = cl | C*cl;
                cl = find(cl);
            else
                cl = opt.W(i).cells;
            end
            if numel(cl) > 1
                cr = cl([end, 1:end-1]);
                I0 = [I0; cl; cr];
                J0 = [J0; cr; cl];
            end
        end
    end
    [order, M, cycleSize] = perform_topo_sort(G, I0, J0, countercurrent);
    cyclesPresent = any(cycleSize > 1);
end

%-------------------------------------------------------------------------%
function [order, M, cycle_size] = perform_topo_sort(G, I0, J0, countercurrent)
    nc = G.cells.num;
    M = sparse(I0, J0, 1, nc, nc);
    processCycles = any(countercurrent);
    if ~processCycles
        order = topological_order(M);
        if isempty(order)
            processCycles = true;
        end
    end
    if processCycles
        [comp, counts] = components(M);
        cycles = counts > 1;

        dispif(mrstVerbose, 'Found %d cycles\n', nnz(cycles));
        [C, ic1, ic2] = unique(comp);

        I1 = ic2(I0);
        J1 = ic2(J0);
        keep = I1 ~= J1;
        I1 = I1(keep);
        J1 = J1(keep);
        m = max(ic2);
        A = sparse(I1, J1, 1, m, m);
        order = topological_order(A);
    else
        dispif(mrstVerbose, 'No cycles found, using topo_sort\n');
        A = M;
        ic2 = ':';
    end
    nn = size(A, 1);
    
    tmp = (1:nn)';
    tmp(order) = 1:nn;
    order = tmp;
    
    order = order(ic2);
    if processCycles
        cycle_size = accumarray(order, 1);
    else
        cycle_size = ones(nc, 1);
    end
end
