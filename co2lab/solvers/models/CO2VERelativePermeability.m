classdef CO2VERelativePermeability < StateFunction

    properties
    end
    
    methods
        function gp = CO2VERelativePermeability(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'s', 'pressure'}, 'state');
            gp.label = 'k_\alpha';
        end

        function kr = evaluateOnDomain(prop, model, state)
            sG = model.getProp(state, 'sg');
            sW = 1 - sG;
            p = model.getProp(state, 'pressure');
            krW = prop.evaluateFluid(model, 'krW', sW, p);
            krG = prop.evaluateFluid(model, 'krG', sG, p);
            kr = {krW, krG};
        end
        
        
    end
    
end
