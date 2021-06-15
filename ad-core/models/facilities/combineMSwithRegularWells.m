function W = combineMSwithRegularWells(W_regular, W_ms)
% Combine regular and MS wells, accounting for missing fields

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

    allflds = unique([fieldnames(W_ms); fieldnames(W_regular)]);
    for i = 1:numel(allflds)
        fld = allflds{i};
        for j = 1:numel(W_regular)
            if ~isfield(W_regular(j), fld)
                W_regular(j).(fld) = [];
            end
        end
        for j = 1:numel(W_ms)
            if ~isfield(W_ms(j), fld)
                W_ms(j).(fld) = [];
            end
        end
    end

    W = [W_regular; W_ms];
    for i = 1:numel(W)
        W(i).isMS = i > numel(W_regular);
    end
end
