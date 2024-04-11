classdef MiscibleWaterComponent < GenericComponent
    
    properties
    end
    
    methods
        function c = MiscibleWaterComponent(name)
            c@GenericComponent(name);
            c = c.functionDependsOn('getComponentDensity', ...
                                    {'ShrinkageFactors', 'SurfaceDensity'}, ...
                                    'PVTPropertyFunctions');
        end
        
        function c = getComponentDensity(component, model, state, varargin)
            c = getComponentDensity@GenericComponent(component, model, state);
            
            % this gives the density of the pure water phase, which in our
            % black-oil like formulation equals the component density of water
            [b, rhoS] = model.getProps(state, 'ShrinkageFactors', 'SurfaceDensity');
            
            phase_ix = model.getPhaseIndex('W');
            
            % component density of water in water phase equals the density of
            % pure water in the water phase 
            c{phase_ix} = rhoS{phase_ix} .* b{phase_ix}; 
        end
        
        function c = getPhaseCompositionSurface(component, model, state, varargin)
            c = getPhaseCompositionSurface@GenericComponent(component, model, state);
            c{model.getPhaseIndex('W')} = 1;
        end
        
        function c = getPhaseComponentFractionInjection(component, model, state, force)
        % Get the volume fraction of the component in each phase (when
        % injecting from outside the domain)
            c = cell(2, 1); % one per phase
            if isfield(force, 'compi')
                comp_i = vertcat(force.compi);
            else
                comp_i = vertcat(force.sat);
            end
            index = model.getPhaseIndex('W');
            ci = comp_i(:, index);
            if any(ci ~= 0)
                c{index} = ci;
            end
            
        end
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
