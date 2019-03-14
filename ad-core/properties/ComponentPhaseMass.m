classdef ComponentPhaseMass < GridProperty
    properties

    end
    
    methods
        function gp = ComponentPhaseMass(varargin)
            gp@GridProperty(varargin{:});
        end
        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            v = cell(ncomp, nph);
            for c = 1:ncomp
                v(c, :) = model.Components{c}.getComponentMass(model, state);
            end
        end
    end
end