classdef ComponentPhaseMassTotalSaturation < ComponentPhaseMass
    properties

    end
    
    methods
        function gp = ComponentPhaseMassTotalSaturation(model, varargin)
            gp@ComponentPhaseMass(model, varargin{:});
            gp = gp.dependsOn('TotalSaturation');
        end
        function v = evaluateOnDomain(prop, model, state)
            v = evaluateOnDomain@ComponentPhaseMass(prop, model, state);
            sT = prop.getEvaluatedDependencies(state, 'TotalSaturation');
            for i = 1:numel(v)
                if ~isempty(v{i})
                    v{i} = v{i}.*sT;
                end
            end
        end
    end
end