function [ws, T] = convertJutulWellSols(wells)
% Convert Jutul wells to MRST wellSols
%
% SYNOPSIS:
%   [ws, T] = convertJutulWellSols(wells)
%
% REQUIRED PARAMETERS:
%   wells - Well output from wells.mat file that Jutul writes.
%
% RETURNS:
%   ws - wellSols that mimick output from MRST's AD solvers
%
%   T  - Reservoir time since simulation start for each of the well sols

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
    T = wells.time;
    names = setdiff(fieldnames(wells), 'time');
    n = numel(T);
    nw = numel(names);
    ws = cell(n, 1);
    
    first = names{1};
    has_water = isfield(wells.(first), 'wrat');
    has_oil = isfield(wells.(first), 'orat');
    has_gas = isfield(wells.(first), 'grat');

    default = struct('name', 'dummy', 'bhp', 0, 'status', true);
    if has_water; default.qWs = 0; end
    if has_oil; default.qOs = 0; end
    if has_gas; default.qGs = 0; end

    ws0 = repmat(default, 1, nw);
    for i = 1:n
        w = ws0;
        for j = 1:nw
            name = names{j};
            W = wells.(name);
            w(j).name = name;
            w(j).bhp = W.bhp(i);
            if has_water
               w(j).qWs = W.wrat(i);
            end
            if has_oil
               w(j).qOs = W.orat(i);
            end
            if has_gas
               w(j).qGs = W.grat(i);
            end
        end
        ws{i} = w;
    end
end
