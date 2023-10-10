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
