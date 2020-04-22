classdef ConcentrationComponent < GenericComponent
    % A component description that assumes that the component is
    % transported in another phase and affect the property of the phase,
    % but it does not change the mass and density of the phase.
    properties
        phaseIndex % Index of phase this component is transported in
    end

    methods
        function c = ConcentrationComponent(name, phase)
            c@GenericComponent(name);
            c.phaseIndex = phase;
            c.isConcentration = true;
            c = c.functionDependsOn('getComponentDensity', 'Density', 'PVTPropertyFunctions');
        end
       
        function c = getPhaseComponentFractionInjection(component, model, state, force)
            % Get the fraction of the component in each phase (when
            % injecting from outside the domain)
            c = cell(model.getNumberOfPhases(), 1);
            if isfield(force, 'compi')
                comp_i = vertcat(force.compi);
            else
                comp_i = vertcat(force.sat);
            end
            index = component.phaseIndex;
            c_inj = component.getInjectionConcentration(force);
            % TODO: the follwing code is using water density explicitly,
            % how to generialize it
            ci = comp_i(:, index) .* c_inj./model.fluid.rhoWS;
            if any(ci ~= 0)
                c{index} = ci;
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
