classdef FixedTotalFlux < StateFunction
    properties(Access = protected)
        fluxfield = 'flux'; % Name of field in state where static flux is stored
    end
    
    methods
        function gp = FixedTotalFlux(model, varargin)
            if nargin > 1 && ischar(varargin{1})
                fld = varargin{1};
                varargin = varargin(2:end);
            else
                fld = 'flux';
            end
            gp@StateFunction(model, varargin{:});
            
            gp.fluxfield = fld;
            gp = gp.dependsOn({fld}, 'state');
        end
        function vT = evaluateOnDomain(prop, model, state)
            vT = sum(state.(prop.fluxfield)(model.operators.internalConn, :), 2);
        end
    end
end