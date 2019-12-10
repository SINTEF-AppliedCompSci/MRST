classdef FixedTotalVelocityDG < FixedTotalFluxDG
    
    properties
    end
    
    methods
        function gp = FixedTotalVelocityDG(model, varargin)
            gp@FixedTotalFluxDG(model, varargin{:});
        end
        function vT = evaluateOnDomain(prop, model, state)
            vT = model.discretization.velocityInterp.faceFlux2cellVelocity(sum(state.flux,2));
            vT = vT(state.cells,:);
        end
    end
end