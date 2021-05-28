classdef SurfactantCapillaryPressure < BlackOilCapillaryPressure
% Scale the water-oil capillary pressure as follows
%     pcow(S_w,c) = pcow(S)*ift(c)/ift(0)
% where ift(c) is the interfacial tension at surfactant concentration c.
    properties
    end
    
    methods
        function prop = SurfactantCapillaryPressure(model, varargin)
            prop@BlackOilCapillaryPressure(model, varargin{:});
            prop = prop.dependsOn('surfactant', 'state');
            assert(model.water && isfield(model.fluid,'ift'))
        end
        
        function pc = evaluateOnDomain(prop, model, state)
            pc = evaluateOnDomain@BlackOilCapillaryPressure(prop, model,state);
            ph = model.getPhaseNames(); iW = find(ph=='W');
            if ~isempty(pc{iW})
                c  = model.getProps(state, 'surfactant');
                pc{iW} = pc{iW}.*model.fluid.ift(c)./model.fluid.ift(0);
            end
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
