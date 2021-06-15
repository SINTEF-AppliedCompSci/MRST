classdef UpwindDiscretization
    % Base class for upwind discretization. The upwind discretization is in
    % general used for the hyperbolic part of equations (e.g. upwinding the
    % mass in an advenction equation)
    methods
        function up = UpwindDiscretization(model)
            
        end
        
        function facevalues = faceUpstream(wrapper, model, state, flag, cellvalues)
            error('Base class not intended for direct usage');
        end
    end
end

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
