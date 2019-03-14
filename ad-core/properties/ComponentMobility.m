classdef ComponentMobility < GridProperty
    properties

    end
    
    methods
        function gp = ComponentMobility(varargin)
            gp@GridProperty(varargin{:});
        end
        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            v = cell(ncomp, nph);
            for i = 1:ncomp
                v(i, :) = model.Components{i}.getComponentMobility(model, state);
            end
        end
    end
end