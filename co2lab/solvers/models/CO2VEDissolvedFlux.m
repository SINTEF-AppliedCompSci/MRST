classdef CO2VEDissolvedFlux < StateFunction & UpwindProperty

    properties (Access = protected)
        upwind_name; % Name of state function where upwind flag comes from
    end

    methods
        function fm = CO2VEDissolvedFlux(model, upstr, upwind_name)
            if nargin < 2
                upstr = UpwindFunctionWrapperDiscretization(model);
            end
            if nargin < 3
                upwind_name = 'PhaseUpwindFlag';
            end
            
            fm@StateFunction(model);
            fm@UpwindProperty(upstr);
            fm.upwind_name = upwind_name;
            fm = fm.dependsOn(upwind_name);
            fm = fm.dependsOn('Mobility', 'FlowPropertyFunctions');
            fm = fm.dependsOn('ShrinkageFactors', 'PVTPropertyFunctions');
            fm = fm.dependsOn('PermeabilityPotentialGradient');
            fm = fm.dependsOn('rs', 'state');
            fm.label = 'v_{rs}';
        end

        function v = evaluateOnDomain(prop, model, state)
            flag = prop.getEvaluatedDependencies(state, prop.upwind_name);
            [mob, b] = prop.getEvaluatedExternals(model, state, 'Mobility', 'ShrinkageFactors');
            
            
            kgrad = prop.getEvaluatedDependencies(state, 'PermeabilityPotentialGradient');
            rs = model.getProp(state, 'rs');

            % @@ Water phase is hard-coded as appearing in position 1 for now
            w_ix = 1;
            
            mob_rs = mob{w_ix} .* rs .* b{w_ix};
            
            v = -1 * prop.faceUpstream(model, state, flag{w_ix}, mob_rs) .* kgrad{w_ix};
        end
        
    end
end
