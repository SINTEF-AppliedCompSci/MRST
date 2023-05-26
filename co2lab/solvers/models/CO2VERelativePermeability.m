classdef CO2VERelativePermeability < BaseRelativePermeability

    properties
    end
    
    methods
        function gp = CO2VERelativePermeability(model, varargin)
            gp@BaseRelativePermeability(model, varargin{:});
            %gp = gp.dependsOn({'s', 'pressure'}, 'state');
            %gp.label = 'k_\alpha';
            gp = gp.dependsOn({'pressure'}, 'state');
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
