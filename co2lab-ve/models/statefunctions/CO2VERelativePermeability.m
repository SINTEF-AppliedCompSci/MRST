classdef CO2VERelativePermeability < BaseRelativePermeability

    properties
    end
    
    methods
        function gp = CO2VERelativePermeability(model, varargin)
            gp@BaseRelativePermeability(model, varargin{:});
            
            if model.hysteresis
                gp = gp.dependsOn({'sGmax', 'pressure'}, 'state');
            end
        end

        function kr = evaluateOnDomain(prop, model, state)
            sG = model.getProp(state, 'sg');
            
            if model.hysteresis
                sGmax = model.getProp(state, 'sGmax');
            else
                sGmax = sG;
            end
            sW = 1 - sG;
            p = model.getProp(state, 'pressure');
            krW = prop.evaluateFluid(model, 'krW', sW, p, 'sGmax', sGmax);
            krG = prop.evaluateFluid(model, 'krG', sG, p, 'sGmax', sGmax);
            kr = {krW, krG};
        end
        
        
    end
    
end
