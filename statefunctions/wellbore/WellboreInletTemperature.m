classdef WellboreInletTemperature < StateFunction

    methods

        function T = evaluateOnDomain(prop, model, state)
            
            T = model.getProp(state, 'bht');
    
        end

    end

end