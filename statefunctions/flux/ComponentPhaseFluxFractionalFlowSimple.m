classdef ComponentPhaseFluxFractionalFlowSimple < StateFunction
    properties
    end
    
    methods
        function cf = ComponentPhaseFluxFractionalFlowSimple(model)
            cf@StateFunction(model);
            cf = cf.dependsOn({'TotalFlux'});
            assert(isfield(model.fluid, 'f_w'), 'Fractional flow must be specified for water');
        end

        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            [vT] = prop.getEvaluatedDependencies(state,'TotalFlux');
            assert(ncomp == nph, 'Immiscible assumption');
            sw  = model.getProp(state,'sw');
            f_w = model.fluid.f_w(sw);
            f_w = model.operators.faceUpstr(vT > 0, f_w);
            
            f_o = 1 - f_w;
            
            v = cell(ncomp, nph);
            v{1, 1} = f_w.*vT;
            v{2, 2} = f_o.*vT;
        end
    end
end
