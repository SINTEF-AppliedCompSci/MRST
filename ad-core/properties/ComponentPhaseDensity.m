classdef ComponentPhaseDensity < GridProperty
    properties

    end
    
    methods
        function gp = ComponentPhaseDensity(model, varargin)
            gp@GridProperty(model, varargin{:});
            
            ncomp = model.getNumberOfComponents();
            deps = cell(ncomp, 1);
            for c = 1:ncomp
                deps{c} = model.Components{c}.dependencies;
            end
            deps = unique(vertcat(deps{:}));
            gp = gp.dependsOn(deps);
        end
        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            v = cell(ncomp, nph);
            for c = 1:ncomp
                v(c, :) = model.Components{c}.getComponentDensity(model, state);
            end
        end
    end
end