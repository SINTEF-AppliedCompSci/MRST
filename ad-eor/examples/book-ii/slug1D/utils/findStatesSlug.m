function [is,ip] = findStatesSlug(cs, cp, slug)
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

i1 = find(cs>.001,1,'first');
i2 = find(cs>.999,1,'first');
i3 = find(cs>1.001,1,'first');
i4 = find(cs>49.99,1,'first');
i5 = find(cs>0.001,1,'last');
is.ix=[1 i1; 1 i2; 1 i3; 1 i4; 1 i5];

j1 = find(cp>.001,1,'first');
j2 = find(cp>2.999,1,'first');
j3 = find(cp>0.001,1,'last');
ip.ix = [1 j1; 1 j2; 1 j3];

switch slug
    case 1
        is.y = [15 30; 15 30; 15 30; 15 30; 5 30];
        ip.y = [15 30; 15 30;  5 30];
    case 2
        is.y = [10 30; 10 30; 10 30; 15 30; 5 30];
        ip.y = [15 30; 15 30; 10 30];
    otherwise
        disp('Not valid value');
end
end
