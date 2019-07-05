classdef ComponentMobilityTotalSaturation < ComponentMobility
    % Class implementing the mobility for a specific component
    properties

    end
    
    methods
        function gp = ComponentMobilityTotalSaturation(model, varargin)
            gp@ComponentMobility(model, varargin{:});
            gp = gp.dependsOn('s', 'state');
        end
        function v = evaluateOnDomain(prop, model, state)
            v = evaluateOnDomain@ComponentMobility(prop, model, state);
            s = model.getProp(state, 's');
            nph = numel(s);
            sT = 0;
            for i = 1:nph
                sT = sT + s{i};
            end
            for ph = 1:nph
                for i = 1:size(v, 1)
                    v{i, ph} = v{i, ph}.*sT;
                end
            end
        end
    end
end