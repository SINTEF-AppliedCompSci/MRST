classdef PhaseUpwindFlagTotalVelocity < StateFunction
    % Upwind flag based on total velocity only (for hybrid upwind)
    properties

    end
    
    methods
        function gp = PhaseUpwindFlagTotalVelocity(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('TotalFlux');
        end

        function flags = evaluateOnDomain(prop, model, state)
            vT = prop.getEvaluatedDependencies(state, 'TotalFlux');
            nph = model.getNumberOfPhases();
            flag = vT > 0;
            flags = cell(1, nph);
            for i = 1:nph
                flags{i} = flag;
            end
        end
    end
end