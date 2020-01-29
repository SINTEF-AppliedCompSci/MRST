classdef FlowPropertyFunctions < StateFunctionGrouping
    % Default grouping for describing a system of flow equations. Contains
    % basic properties like mobility, density, component total mass etc.
    properties
        RelativePermeability % Phase relative permeabilities (1 x nphase)
        CapillaryPressure % Capillary pressures (if present) (1 x nphase, with empty entries for phases with equal pressure to reference)
        Mobility % Phase mobilities (1 x nphase)
        % Component properties
        ComponentTotalMass % Total component mass (ncomp x 1)
        ComponentPhaseMass % Component mass in each phase (ncomp x nphase)
        ComponentMobility % Mobility of each component in each phase (mass, not volume) (ncomp x nphase)
        ComponentPhaseDensity % Mass-density per volume in each phase (ncomp x nphase)
    end

    methods
        function props = FlowPropertyFunctions(model)
            props@StateFunctionGrouping('FlowProps');
            sat = props.getRegionSaturation(model);
            % Saturation properties
            props = props.setStateFunction('CapillaryPressure', BlackOilCapillaryPressure(model, sat));
            props = props.setStateFunction('RelativePermeability', BaseRelativePermeability(model, sat));
            props = props.setStateFunction('Mobility', Mobility(model, sat));

            % Components
            props = props.setStateFunction('ComponentPhaseMass', ComponentPhaseMass(model));
            props = props.setStateFunction('ComponentTotalMass', ComponentTotalMass(model));
            props = props.setStateFunction('ComponentMobility', ComponentMobility(model));
            props = props.setStateFunction('ComponentPhaseDensity', ComponentPhaseDensity(model));
        end
        
        function sat = getRegionSaturation(props, model)
            r = model.rock;
            sat = ones(model.G.cells.num, 1);
            if isfield(r, 'regions')
                if isfield(r.regions, 'saturation')
                    sat = r.regions.saturation;
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
