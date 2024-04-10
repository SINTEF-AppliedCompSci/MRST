classdef WellboreBottomholeTemperature < StateFunction
% State function for bottom-hole temperature in WellboreModel.m
    
    methods
        %-----------------------------------------------------------------%
        function cf = WellboreBottomholeTemperature(model)

            cf@StateFunction(model);
            cf = cf.dependsOn({'T'}, 'state');
            cf.label = 'T_{bh}';

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function bht = evaluateOnDomain(prop, model, state)
            
            % Bottom hole temperature is defined to be the pressure in the
            % inlet segment.
            T   = model.getProps(state, 'T');
            iic = model.getInletSegments();
            bht = T(iic);
            
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