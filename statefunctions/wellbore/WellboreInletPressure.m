classdef WellboreInletPressure < StateFunction

    methods

        function p = evaluateOnDomain(prop, model, state)
            
            p = model.getProp(state, 'bhp');
    
        end

    end

end