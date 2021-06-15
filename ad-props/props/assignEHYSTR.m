function f = assignEHYSTR(f, ehystr, reg)
%
% Assigns flags for hysteresis computation, as well as input options for
% hysteresis model according to ECLIPSE EHYSTR keyword. 

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

% flags for hysteresis
if strcmp(ehystr{5}, 'KR') || strcmp(ehystr{5}, 'BOTH')
    f.krHyst = 1;
elseif strcmp(ehystr{5}, 'PC') || strcmp(ehystr{5}, 'BOTH')
    f.pcHyst = 1;
end

% input options
f.ehystr = ehystr;

end
