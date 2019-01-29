classdef FacilityFluxDiscretization < PropertyFunctions
    properties
        PhaseFlux
        ComponentTotalFlux
        ComponentPhaseFlux
        PerforationPressureGradient
        WellIndex
    end
    
    methods
        function props = FacilityFluxDiscretization(model)
            props.structName = 'FacilityFluxProps';
            
            props.PhaseFlux = WellPhaseFlux(model);
            props.ComponentTotalFlux = ComponentTotalFlux(model);
            props.ComponentPhaseFlux = WellComponentPhaseFlux(model);
            props.PerforationPressureGradient = PerforationPressureGradient(model);
            props.WellIndex = WellIndex(model);
        end
    end
end