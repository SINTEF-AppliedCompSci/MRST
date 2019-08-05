classdef ComponentPhaseDensityTotalSaturation < ComponentPhaseDensity
    properties

    end
    
    methods
        function gp = ComponentPhaseDensityTotalSaturation(model, varargin)
            gp@ComponentPhaseDensity(model, varargin{:});
            gp = gp.dependsOn('TotalSaturation');
        end
        function v = evaluateOnDomain(prop, model, state)
            v = evaluateOnDomain@ComponentPhaseDensity(prop, model, state);
            sT = prop.getEvaluatedDependencies(state, 'TotalSaturation');
            for i = 1:numel(v)
                if ~isempty(v{i})
                    v{i} = v{i}.*sT;
                end
            end
        end
    end
end