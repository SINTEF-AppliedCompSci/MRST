classdef PerforationComponentPhaseDensity < StateFunction
    % Component density to used for each well connection
    properties

    end
    
    methods
        function gp = PerforationComponentPhaseDensity(varargin)
            gp = gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'ComponentPhaseDensity'}, 'FlowPropertyFunctions');
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp.label = '\rho_{wc}';
        end
        
        function rhoc = evaluateOnDomain(prop, model, state)
            rhoc = model.ReservoirModel.getProps(state, 'ComponentPhaseDensity');
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            for i = 1:numel(rhoc)
                if ~isempty(rhoc{i})
                    rhoc{i} = rhoc{i}(map.cells);
                end
            end
        end
    end
end
