classdef Temperature < StateFunction
   
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = Temperature(model, varargin)
            gp@StateFunction(model, varargin{:});
            switch model.thermalFormulation
                case 'enthalpy'
                    % Enthalpy is the primary variable
                    gp = gp.dependsOn({'PhasePressures', 'Density', 'Enthalpy'});
                case 'temperature'
                    gp = gp.dependsOn({'temperature'}, 'state');
                    % Temperature is the primary variable - we can get it
                    % directly from state
            end
            gp.label = 'T';
        end
        
        %-----------------------------------------------------------------%
        function T = evaluateOnDomain(prop, model, state)
            switch model.thermalFormulation
                case 'enthalpy'
                    error('Not implemented yet')
                case 'temperature'
                    T = state.T;
            end
        end
    end
    
end