function CG = basisToCoarseGrid(G, I)
% Create a coarse grid from a set of basis functions so that each coarse
% block corresponds to the places where that basis is the largest value.

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

    for i = 1:size(I, 2)
        [~, biggest] = max(I(:, i));
        I(biggest, i) = inf;
    end
    [~, p] = max(I, [], 2);
    
    CG = generateCoarseGrid(G, p);
    CG = coarsenGeometry(CG);
end
