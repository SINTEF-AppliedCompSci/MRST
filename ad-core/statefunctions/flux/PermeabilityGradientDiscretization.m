classdef PermeabilityGradientDiscretization
    % Class defining the discretization of the permeability
    methods
        function facevalues = getPermeabilityGradient(tpfa, model, state, facepotentialdiff)
            % Should return the discretized version of K * dp on an
            % interior face for a potential difference dp (also defined on
            % internal faces).
        end
    end
end