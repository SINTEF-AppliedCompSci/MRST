classdef BrineComponent < GenericComponent
%Class for components in brine
    
    properties
        componentIndex
        molecularDiffusivity
        surfacePhaseMassFractions = 1;
    end
    
    methods
        %-----------------------------------------------------------------%
        function component = BrineComponent(name, molarMass, diffusivity, index)
            component                      = component@GenericComponent(name);
            component.molarMass            = molarMass;
            component.molecularDiffusivity = diffusivity;
            component.componentIndex       = index;
        end
        
        %-----------------------------------------------------------------%
        function componentDensity = getComponentDensity(component, model, state, extra)
            % Density of component in each phase (mass per unit of volume)
            [rho, X] = model.getProps(state, 'Density', 'ComponentPhaseMassFractions');
            ix = strcmp(model.getComponentNames, component.name);
            rho{1} = rho{1}.*X{ix};
            componentDensity = rho;
        end
        
        %-----------------------------------------------------------------%
        function c = getPhaseComponentFractionInjection(component, model, state, force)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            comp_i = vertcat(force.components);
            comp_i = model.getMassFraction(comp_i);
            index = component.componentIndex;
            ci = comp_i(:, index);
            if any(ci ~= 0)
                for i = 1:nph
                    c{i} = ci;
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function c = getPhaseCompositionSurface(component, model, state, pressure, temperature)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            for i = 1:nph
                c{i} = component.surfacePhaseMassFractions(i);
            end
        end
    end
    
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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