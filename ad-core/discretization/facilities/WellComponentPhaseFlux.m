classdef WellComponentPhaseFlux < GridProperty
    properties

    end
    
    methods

        function cflux = evaluateOnDomain(prop, facility, state)
            model = facility.ReservoirModel;
            q = facility.getProp(state, 'PhaseFlux');
            rhoc = model.getProp(state, 'ComponentPhaseDensity');
        end
    end
end