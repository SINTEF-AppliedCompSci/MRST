classdef ComponentMobility < StateFunction & ComponentProperty
    % Class implementing the mobility for a specific component
    properties

    end
    
    methods
        function gp = ComponentMobility(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp@ComponentProperty(model);
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