classdef CO2VEDissolvedFlux < StateFunction 

    properties (Access = protected)
    end

    methods
        function fm = CO2VEDissolvedFlux(model)
            fm@StateFunction(model);
            fm = fm.dependsOn('ComponentPhaseFlux');
            fm.label = 'v_{rs}';
        end

        function v = evaluateOnDomain(prop, model, state)
            c_phase_fluxes = prop.getEvaluatedDependencies(state, 'ComponentPhaseFlux');
            
            water_phase_ix = model.getPhaseIndex('W');
            co2_phase_ix = model.getPhaseIndex('G');
            
            water_component_ix = 1; % @@ hard-coded for now
            co2_component_ix = 2;   % @@
            
            v = c_phase_fluxes{co2_component_ix, water_phase_ix};
        end
    end
end

