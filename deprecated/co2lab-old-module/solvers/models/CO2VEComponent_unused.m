classdef CO2VEComponent < ImmiscibleComponent

    methods
        
        function c = CO2VEComponent()
            name = 'CO2';
            phase = 2;
            c@ImmiscibleComponent(name, phase)
        end
        
        function mass = getComponentMass(component, model, state, varargin)
        % @@@ IMPLEMENT ME
        end
        
    end
end
