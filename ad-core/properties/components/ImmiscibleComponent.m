classdef ImmiscibleComponent < ComponentImplementation
    properties
        phaseIndex % Index of phase this component belongs to
    end
    
    methods
        function c = ImmiscibleComponent(name, phase)
            c@ComponentImplementation(name);
            c.phaseIndex = phase;
        end
        
        function c = getComponentDensity(component, model, state, varargin)
            c = getComponentDensity@ComponentImplementation(component, model, state, varargin{:});
            rho = model.getProp(state, 'Density');
            c{component.phaseIndex} = rho{component.phaseIndex};
        end
        
        function c = getPhaseComposition(component, model, state, varargin)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            c{component.phaseIndex} = 1;
        end
        
        function c = getPhaseCompositionWell(component, model, state, W)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            comp_i = vertcat(W.compi);
            for i = 1:nph
                c{i} = comp_i(:, i);
            end
        end
    end
end