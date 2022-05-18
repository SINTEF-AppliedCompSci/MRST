classdef RockMass < StateFunction

    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = RockMass(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'RockDensity', 'PoreVolume'});
            gp.label = 'M_R';
        end
        
        %-----------------------------------------------------------------%
        function massR = evaluateOnDomain(prop,model, state)
            [rhoR, pv] = prop.getEvaluatedDependencies(state, 'RockDensity', ...
                                                              'PoreVolume' );
            v     = model.operators.vol - pv;
            massR = rhoR.*v;
        end       
    end
    
end