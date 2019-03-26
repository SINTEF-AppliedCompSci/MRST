classdef ComponentPhaseMass < GridProperty & ComponentProperty
    properties

    end
    
    methods
        function gp = ComponentPhaseMass(model, varargin)
            gp@GridProperty(model, varargin{:});
            gp@ComponentProperty(model);
            gp = gp.dependsOn('Components', 'model');
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