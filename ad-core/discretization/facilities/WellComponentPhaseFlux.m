classdef WellComponentPhaseFlux < GridProperty
    properties

    end
    
    methods

        function cflux = evaluateOnDomain(prop, facility, state)
            model = facility.ReservoirModel;
            q = facility.getProp(state, 'PhaseFlux');
        end
    end
end