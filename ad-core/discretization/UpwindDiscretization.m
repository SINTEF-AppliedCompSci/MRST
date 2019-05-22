classdef UpwindDiscretization
    % Base class for upwind discretization. The upwind discretization is in
    % general used for the hyperbolic part of equations (e.g. upwinding the
    % mass in an advenction equation)
    methods
        function up = UpwindDiscretization()
            
        end
        
        function facevalues = faceUpstream(wrapper, state, flag, cellvalues)
            error('Base class not intended for direct usage');
        end
    end
end