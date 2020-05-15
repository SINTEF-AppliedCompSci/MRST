classdef WellConCompPhaseDensity < StateFunction
    % Component density to used for each well connection
    properties

    end
    
    methods
        function gp = WellConCompPhaseDensity(varargin)
            gp = gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'ComponentPhaseDensity'}, 'FlowPropertyFunctions');
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp.label = '\rho_{wc}';
        end
        
        
        function rhoc = evaluateOnDomain(prop, model, state)
            componentPhaseDensity = model.ReservoirModel.getProps(state, 'ComponentPhaseDensity');
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            function v = ex(x, cells)
                v = [];
                if ~isempty(x) 
                    v=x(cells); 
                end
            end
            f = @(x) ex(x, map.cells);
            rhoc = cellfun(f, componentPhaseDensity, 'UniformOutput', false); 
        end
    end
    
end
