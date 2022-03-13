function subtraps = getMigrationTree(G, A, trap, depth)
% Recursively traverse and find the full migration tree
    subtraps = find(A(trap, :));
    tmp = [];
    for i = 1:numel(subtraps);
        trp = getMigrationTree(G, A, subtraps(i), depth + 1);
        tmp = [tmp trp]; %#ok
    end
    subtraps = [reshape(subtraps, 1, []) reshape(tmp, 1, [])];
end

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
