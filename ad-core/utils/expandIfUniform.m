function v = expandIfUniform(v)
    % Utility which reverses "value" compaction. If given a matrix (logical
    % or numerical) as input, it will expand it to a cell array of vectors
    % such that value(expandIfUniform(x)) is equal to x.

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

    if (isnumeric(v) || islogical(v)) && size(v, 2) > 1
        n = size(v, 2);
        out = cell(1, n);
        for i = 1:n
            out{i} = v(:, i);
        end
        v = out;
    end
end
