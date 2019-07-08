classdef PhaseInterfacePressureDifferences < StateFunction
    properties

    end
    
    methods
        function gp = PhaseInterfacePressureDifferences(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('CapillaryPressure', 'FlowPropertyFunctions');
            gp = gp.dependsOn({'PhasePotentialDifference'});
        end

        function G = evaluateOnDomain(prop, model, state)
            pot = prop.getEvaluatedDependencies(state, ...
                            'PhasePotentialDifference');
            pc = model.getProp(state, 'CapillaryPressure');
            G = pot;
            for i = 1:numel(G)
                G{i} = -pot{i};
                if ~isempty(pc{i})
                    G{i} = G{i} + model.operators.Grad(pc{i});
                end
            end
        end
    end
end