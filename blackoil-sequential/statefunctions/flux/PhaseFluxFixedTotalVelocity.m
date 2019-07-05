classdef PhaseFluxFixedTotalVelocity < StateFunction
    properties
    end
    
    methods
        function gp = PhaseFluxFixedTotalVelocity(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('Mobility', 'FlowPropertyFunctions');
            gp = gp.dependsOn({'FractionalFlow', 'PhaseInterfacePressureDifferences'});
            gp = gp.dependsOn({'flux'}, 'state');
        end
        function f = evaluateOnDomain(prop, model, state)
            [f, G] = prop.getEvaluatedDependencies(state, 'FractionalFlow', 'PhaseInterfacePressureDifferences');
            [mob, flux] = model.getProps(state, 'Mobility', 'flux');
            vT = sum(flux(model.operators.internalConn), 2);
        end
    end
end