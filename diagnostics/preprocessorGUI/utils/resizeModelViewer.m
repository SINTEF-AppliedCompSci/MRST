function resizeModelViewer(src, event)
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

fig = src;
c = fig.Children;
histix = find(strcmp('histogram', get(c,'Tag')));
cbarix = find(strcmp('Colorbar', get(c,'Tag')));

for k = 1:numel(histix)
    c(histix(k)).Position(1) = sum(c(cbarix(k)).Position([1 3])) + 0.002;
end
end
