classdef PhasePressures < AutoDiffFunction
    properties
    end
    
    methods
        function gp = PhasePressures(varargin)
            gp@AutoDiffFunction(varargin{:});
            gp = gp.dependsOn({'CapillaryPressure'});
            gp = gp.dependsOn({'pressure'}, 'state');
        end
        
        function p_phase = evaluateOnDomain(prop, model, state)
            p = model.getProps(state, 'Pressure');
            pc = prop.getEvaluatedDependencies(state, 'CapillaryPressure');
            nph = numel(pc);
            p_phase = cell(1, nph);
            for i = 1:nph
                if isempty(pc{i})
                    p_phase{i} = p;
                else
                    p_phase{i} = p + pc{i};
                end
            end
        end
    end
end