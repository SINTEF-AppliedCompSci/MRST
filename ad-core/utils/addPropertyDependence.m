function prop = addPropertyDependence(prop, name, grouping)
    % Document dependencies and external dependencies

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

    if iscell(name)
        name = reshape(name, [], 1);
    end
    if nargin < 3 || isempty(grouping)
        if isstruct(name)
            prop.externals = [prop.externals; name];
        else
            prop.dependencies = [prop.dependencies; name];
        end
    else
        s = struct('name', name, 'grouping', grouping);
        prop.externals = [prop.externals; s];
    end
end
