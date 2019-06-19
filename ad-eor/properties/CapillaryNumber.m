classdef CapillaryNumber < StateFunction
    
    properties
    end
    
    methods
        function prop = CapillaryNumber(varargin)
            prop@StateFunction(varargin{:});
        end
        function Nc = evaluateOnDomain(prop, model, state)
        % For this grid property, it is expected that it has been set before.
            assert(structPropEvaluated(state, 'CapillaryNumber'));
            Nc = state.CapillaryNumber;
        end
    end
end
