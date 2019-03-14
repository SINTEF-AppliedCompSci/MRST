classdef ComponentPhaseDensity < GridProperty
    properties

    end
    
    methods
        function gp = ComponentPhaseDensity(varargin)
            gp@GridProperty(varargin{:});
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