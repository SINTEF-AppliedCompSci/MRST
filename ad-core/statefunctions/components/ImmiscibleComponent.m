classdef ImmiscibleComponent < GenericComponent
    % Specialized interface for an immiscible component
    %
    % The component description assumes that the component is immiscible,
    % i.e. it only exists in one phase that is made up entirely of that
    % specific component.
    properties
        phaseIndex % Index of phase this component belongs to
    end
    
    methods
        function c = ImmiscibleComponent(name, phase)
            c@GenericComponent(name);
            c.phaseIndex = phase;
            c = c.functionDependsOn('getComponentDensity', 'Density', 'PVTPropertyFunctions');
        end
        
        function c = getComponentDensity(component, model, state)
            rho = model.PVTPropertyFunctions.get(model, state, 'Density');
            c = getComponentDensity@GenericComponent(component, model, state);
            if numel(c) == 1 && ~iscell(rho)
                c{1} = rho;
            else
                c{component.phaseIndex} = rho{component.phaseIndex};
            end
        end
        
        function c = getPhaseComposition(component, model, state, varargin)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            c{component.phaseIndex} = 1;
        end
        
        function c = getPhaseCompositionSurface(component, model, state, pressure, temperature)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            c{component.phaseIndex} = 1;
        end
        
        function c = getPhaseComponentFractionInjection(component, model, state, force)
            % Get the volume fraction of the component in each phase (when
            % injecting from outside the domain)
            c = cell(model.getNumberOfPhases(), 1);
            if isfield(force, 'compi')
                comp_i = vertcat(force.compi);
            else
                comp_i = vertcat(force.sat);
            end
            index = component.phaseIndex;
            ci = comp_i(:, index);
            if any(ci ~= 0)
                c{index} = ci;
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
