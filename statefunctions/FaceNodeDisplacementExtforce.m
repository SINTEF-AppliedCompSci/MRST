classdef FaceNodeDisplacementExtforce < StateFunction

    methods
        function gp = FaceNodeDisplacementExtforce(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'displacement', 'lambdamech', 'pressure', 'extforce'}, 'state');
        end
        
        function fndisp = evaluateOnDomain(prop, model, state)
            
            fndispop = model.operators.facenodedispop;
            [u, p, lm, ef] = model.getProps(state, 'displacement', 'pressure', 'lambdamech', 'extforce');
            fndisp = fndispop(u, p, lm, ef);
            
        end
    end
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
