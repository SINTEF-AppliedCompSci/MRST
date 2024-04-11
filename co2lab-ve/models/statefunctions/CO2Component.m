classdef CO2Component < GenericComponent

    properties
    end
    
    methods
        function c = CO2Component(name, disgas)
            c@GenericComponent(name);
            c = c.functionDependsOn('getComponentDensity', ...
                                    {'ShrinkageFactors', 'SurfaceDensity'}, ...
                                    'PVTPropertyFunctions');
            if disgas
                c = c.functionDependsOn('getComponentDensity', 'rs', 'state');
            end
        end
        
        function c = getPhaseCompositionSurface(component, model, state, varargin)
            c = getPhaseCompositionSurface@GenericComponent(component, model, state);
            c{model.getPhaseIndex('G')} = 1;
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
            index = model.getPhaseIndex('G');
            ci = comp_i(:, index);
            if any(ci ~= 0)
                c{index} = ci;
            end
        end
        
        function c = getComponentDensity(component, model, state, varargin)
            c = getComponentDensity@GenericComponent(component, model, state, varargin{:});
            
            [b, rhoS] = model.getProps(state, 'ShrinkageFactors', ...
                                              'SurfaceDensity');
            co2_phase_ix = model.getPhaseIndex('G');
            
            % component density of CO2 in CO2 phase 
            c{co2_phase_ix} = rhoS{co2_phase_ix} .* b{co2_phase_ix};
            
            if model.disgas
                
                % component density of CO2 in brine phase
                w_phase_ix = model.getPhaseIndex('W');
                
                rs = model.getProp(state, 'rs');
                rhoGS = rhoS{co2_phase_ix};
                bW = b{w_phase_ix};
                
                c{w_phase_ix} = rs .* rhoGS .* bW;
            end
        end
    end
end
