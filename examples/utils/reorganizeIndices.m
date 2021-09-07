function [lumping,subset] = reorganizeIndices(indices)
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

A = cellfun(@(C)numel(C),indices);

sizeAllIndices = sum(A);

lumping = zeros(sizeAllIndices,1);
subset  = zeros(sizeAllIndices,1);

reel = 1;
for i=1:numel(indices)
    n = numel(indices{i});
    subset(reel:n+reel-1) = indices{i};
    lumping(reel:n+reel-1)= i;
    reel = reel+ n;
end
