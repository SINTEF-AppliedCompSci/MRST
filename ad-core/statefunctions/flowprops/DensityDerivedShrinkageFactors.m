classdef DensityDerivedShrinkageFactors < StateFunction
    properties
    end
    
    methods
        function gp = DensityDerivedShrinkageFactors(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'Density'});
        end

        function b = evaluateOnDomain(prop, model, state)
            b = prop.getEvaluatedDependencies(state, 'Density');
            rhoS = model.getSurfaceDensities();
            for i = 1:numel(b)
                b{i} = b{i}./rhoS(i);
            end
        end
    end
end