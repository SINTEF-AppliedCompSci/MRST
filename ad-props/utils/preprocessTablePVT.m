function tab = preprocessTablePVT(tab)
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

    nk = numel(tab.key);
    exp = cell(1, nk);
    for tn = 1:nk
        % Precompute each table from the compressed format
        lns = tab.pos(tn):(tab.pos(tn+1)-1);
        T = tab.data(lns, :);
        if size(T, 1) > 1 && T(1, 1) > T(2, 1)
            % If the data is decreasing rather than increasing, we reverse
            % the table.
            T = T(end:-1:1, :);
        end
        exp{tn} = struct('x', T(:, 1), 'F', T(:, 2:end));
    end
    tab.expanded = exp;
    T_sat = tab.data(tab.pos(1:end-1), :);
    tab.sat = struct('x', T_sat(:, 1), 'F', T_sat(:, 2:end));
end
