classdef FixedTotalVelocityDG < FixedTotalFluxDG
    
    properties
    end
    
    methods
        function gp = FixedTotalVelocityDG(model, varargin)
            gp@FixedTotalFluxDG(model, varargin{:});
        end
        function vT = evaluateOnDomain(prop, model, state)
            vT = faceFlux2cellVelocity(model.G, sum(state.(prop.fluxfield),2));
            vT = mat2cell(vT, model.G.cells.num, ones(1,model.G.griddim));
            vT = SpatialVector(vT{:});
            vT = vT(state.cells,:);
        end
    end
end