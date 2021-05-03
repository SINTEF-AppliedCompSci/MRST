classdef GenericComponent < StateFunctionDependent
    % The base class that implements all interfaces for a given component.
    % The helper class provides a number of utility functions that
    % determine the mobility, total mass and composition at different
    % conditions.
    properties
        name
        externals = [];
        molarMass = 1;
    end
    
    methods
        function c = GenericComponent(name)
            c.name = name;
            % Document dependencies internal to grouping
            c = c.functionDependsOn('getComponentMobility', {'Mobility'}, 'FlowPropertyFunctions');
            c = c.functionDependsOn('getComponentMobility', {'PoreVolume', 'Density'}, 'PVTPropertyFunctions');
            % State dependencies
            c = c.functionDependsOn('getComponentMass', {'PoreVolume', 'Density'}, 'PVTPropertyFunctions');
            c = c.functionDependsOn('getComponentMass', 's', 'state');
        end
        
        function c = getPhaseComposition(component, model, state)
            % Mass fraction of component in each phase. Default implementation
            % goes via the component density and total density and assumes
            % these are internally consistent.
            c = component.getComponentDensity(model, state);
            rho = model.getProp(state, 'Density');
            if ~iscell(rho)
                rho = {rho};
            end
            for ph = 1:numel(c)
                if ~isempty(c{ph})
                    c{ph} = c{ph}./rho{ph};
                end
            end
        end
        
        function mass = getComponentMass(component, model, state)
            % Mass of component in each phase
            % saturation * pore-volume * component mass density
            s = model.getProp(state, 's');
            s = expandMatrixToCell(s);
            pv = model.PVTPropertyFunctions.get(model, state, 'PoreVolume');
            mass = component.getComponentDensity(model, state);
            % Iterate over phases and weight by pore-volume and saturation
            for i = 1:numel(mass)
                if ~isempty(mass{i})
                    mass{i} = s{i}.*pv.*mass{i};
                end
            end
        end
        
        function cmob = getComponentMobility(component, model, state, varargin)
            % The amount of mobile mass in the cell. For most models, this
            % is just the mobility multiplied with the component density in
            % the cell.
            mob = model.FlowPropertyFunctions.get(model, state, 'Mobility', true);
            density = component.getComponentDensity(model, state, varargin{:});
            nphase = numel(density);
            cmob = cell(1, nphase);
            for i = 1:nphase
                if ~isempty(density{i})
                    cmob{i} = mob{i}.*density{i};
                end
            end
        end
        
        function c = getComponentDensity(component, model, state, extra)
            % Density of component in each phase (mass per unit of volume)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
        end

        function c = getPurePhaseDensitySurface(component, model, state, pressure, temperature)
            % Density of component in each phase at surface conditions
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
        end
        
        function c = getPhaseCompositionSurface(component, model, state, pressure, temperature)
            % Mass fraction of the component in each phase at surface
            % conditions (specified with pressure and temperature)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
        end
        
        function c = getPhaseComponentFractionInjection(component, model, state, force)
            % Get the fraction of the component in each phase (when
            % injecting from outside the domain)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
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
