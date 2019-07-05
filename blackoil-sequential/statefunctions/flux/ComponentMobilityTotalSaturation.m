classdef ComponentMobilityTotalSaturation < ComponentMobility
    % Class implementing the mobility for a specific component
    properties

    end
    
    methods
        function gp = ComponentMobilityTotalSaturation(model, varargin)
            gp@ComponentMobility(model, varargin{:});
            gp = gp.dependsOn('TotalSaturation');
        end
        function v = evaluateOnDomain(prop, model, state)
            v = evaluateOnDomain@ComponentMobility(prop, model, state);
            sT = prop.getEvaluatedDependencies(state, 'TotalSaturation');
            for i = 1:numel(v)
                if ~isempty(v{i})
                    v{i} = v{i}.*sT;
                end
            end
        end
    end
end