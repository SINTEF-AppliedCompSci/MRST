function T = TransNTPFA(G, u, OSflux)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    dispif(mrstVerbose, 'TransNTPFA... ');
    timer = tic;

    T = cell(2, 1);

    internal = 1:G.faces.num;
    internal(~all(G.faces.neighbors ~= 0, 2)) = [];

    tii = cell(2, 1);
    tsp = cell(2, 1);
    r = cell(2, 1);
    mu = cell(2, 1);

    for j = 1:2
        tend = zeros(G.faces.num, 1);
        tii{j} = zeros(G.faces.num, 2);
        nij = zeros(G.faces.num, 1);

        % Sweep for sparsity set up
        for i = internal
            t = OSflux{i, j};            
            nij(i) = max(0, size(t, 1) - 3);
            tii{j}(i, 1:2) = [t(1, 2), t(2, 2)];
        end

        % Set up sparsity pattern and values
        ii = zeros(sum(nij), 1);
        jj = zeros(sum(nij), 1);
        vv = zeros(sum(nij), 1);
        s = [1; cumsum(nij) + 1];

        % Don't loop over zero rows
        nzrows = 1:G.faces.num;
        nzrows(nij == 0) = [];

        for i = nzrows
            t = OSflux{i, j};
            idx = s(i):(s(i+1) - 1);
            ii(idx) = i;
            jj(idx) = t(3:end-1, 1);
            vv(idx) = t(3:end-1, 2);
            tend(i) = t(end, 2);
        end

        tsp{j} = sparse(ii, jj, vv, G.faces.num, G.cells.num);

        r{j} = tsp{j} * u + tend;
        mu{j} = 0 * r{j} + 0.5;
        T{j} = 0 * r{j};
    end

    epstol = 1e-12 * max(full(max(tsp{1}, [], 2)), full(max(tsp{2}, [], 2)));
    for j = 1:2
        ir = abs(r{j}) <= epstol;
        r{j}(ir) = 0;
    end
    jj = abs(r{1} + r{2}) > epstol;
    mu{1}(jj) = r{2}(jj) ./ (r{1}(jj) + r{2}(jj));
    mu{2}(jj) = ones(sum(jj), 1) - mu{1}(jj);
    assert(all(mu{1} >= 0.0), ['min(mu{1})=', num2str(minmax(mu{1}))])
    assert(all(mu{2} >= 0.0), ['min(mu{2})=', num2str(minmax(mu{2}))])
    assert(all(mu{1} <= 1.0), ['max(mu{1})=', num2str(minmax(mu{1}))])
    assert(all(mu{2} <= 1.0), ['max(mu{2})=', num2str(minmax(mu{2}))])

    T{1}(internal) = mu{1}(internal) .* tii{1}(internal, 1) + mu{2}(internal) .* tii{2}(internal, 2);
    T{2}(internal) = mu{1}(internal) .* tii{1}(internal, 2) + mu{2}(internal) .* tii{2}(internal, 1);

    dispif(mrstVerbose, 'done in %1.2f s\n', toc(timer));
end

function y = minmax(x)

    if isa(x, 'ADI')
        y = [min(x.val), max(x.val)];
    else
        y = [min(x), max(x)];
    end

end
