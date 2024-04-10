classdef WellboreBottomholePressure < StateFunction
% State function for bottom-hole pressure in WellboreModel.

    methods
        %-----------------------------------------------------------------%
        function cf = WellboreBottomholePressure(model)

            cf@StateFunction(model);
            cf = cf.dependsOn({'pressure'}, 'state');
            cf.label = 'p_{bh}';

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function bhp = evaluateOnDomain(prop, model, state)
        
            % Bottom hole pressure is defined to be the pressure in the
            % inlet segment.
            p   = model.getProps(state, 'pressure');
            iic = model.getInletSegments();
            bhp = p(iic);
            
        end
        %-----------------------------------------------------------------%
        
    end
    
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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