classdef FixedTotalFlux < StateFunction
    properties(Access = protected)
        fluxfield = 'flux'; % Name of field in state where static flux is stored
    end
    
    methods
        function gp = FixedTotalFlux(model, varargin)
            if nargin > 1 && ischar(varargin{1})
                fld = varargin{1};
                varargin = varargin(2:end);
            else
                fld = 'flux';
            end
            gp@StateFunction(model, varargin{:});
            
            gp.fluxfield = fld;
            gp = gp.dependsOn({fld}, 'state');
        end
        function vT = evaluateOnDomain(prop, model, state)
            vT = sum(state.(prop.fluxfield)(model.operators.internalConn, :), 2);
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
