classdef CO2VEBlackOilDensity < StateFunction

    properties
    end
    
    methods
        function gp = CO2VEBlackOilDensity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'ShrinkageFactors', 'SurfaceDensity'});
            if model.disgas
                gp = gp.dependsOn({'rs'}, 'state');
            end
            gp.label = '\rho_\alpha';
            gp.outputRange = [0, inf];
        end
        
        function rho = evaluateOnDomain(prop, model, state)
            [b, rhoS] = prop.getEvaluatedDependencies(state, 'ShrinkageFactors', ...
                                                             'SurfaceDensity');
            nph = size(rhoS, 2);
            rho = cell(1, nph);

            % First, include the pure phase densities
            for i = 1:nph
                rho{i} = rhoS{i}.*b{i};
            end
            
            % Include dissolved CO2 density in brine if applicable
            if model.disgas
                [wix, gix] = model.getPhaseIndex('W', 'G');
                rs = model.getProp(state, 'rs');
                rho{wix} = rho{wix} + rs.*b{wix}.*rhoS{gix};
            end
        end
    end
end
