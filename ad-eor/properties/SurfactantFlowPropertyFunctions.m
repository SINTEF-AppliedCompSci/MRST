classdef SurfactantFlowPropertyFunctions < FlowPropertyFunctions
    properties
        CapillaryNumber
        SurfactantAdsorption
    end
    
    methods
        function props = SurfactantFlowPropertyFunctions(model)
            props = props@FlowPropertyFunctions(model);
            satreg  = props.getRegionSaturation(model);
            surfreg = props.getRegionSurfactant(model);
            
            props.RelativePermeability = SurfactantRelativePermeability(model, ...
                                                              satreg, surfreg);
            props.CapillaryPressure    = SurfactantCapillaryPressure(model, satreg);
            props.CapillaryNumber      = CapillaryNumber(model);
            props.SurfactantAdsorption = SurfactantAdsorption(model);
            props.Viscosity = BlackOilSurfactantViscosity(model, satreg);
            props.PhasePressures = SurfactantPhasePressures(model, satreg);
        end
        
        function sat = getRegionSurfactant(props, model)
            r = model.rock;
            sat = ones(model.G.cells.num, 1);
            if isfield(r, 'regions')
                if isfield(r.regions, 'surfactant')
                    sat = r.regions.surfactant;
                end
            end
        end
        
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
