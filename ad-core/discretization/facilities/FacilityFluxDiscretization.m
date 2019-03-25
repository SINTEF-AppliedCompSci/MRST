classdef FacilityFluxDiscretization < PropertyFunctions
    properties
        PhaseFlux
        ComponentTotalFlux
        ComponentPhaseFlux
        PerforationPressureGradient
        WellIndex
        FacilityWellMapping
    end
    
    methods
        function props = FacilityFluxDiscretization(model)
            props.structName = 'FacilityFluxProps';
            
            backend = model.AutoDiffBackend;
            props.PhaseFlux = WellPhaseFlux(backend);
            props.ComponentTotalFlux = ComponentTotalFlux(backend);
            props.ComponentPhaseFlux = WellComponentPhaseFlux(backend);
            props.PerforationPressureGradient = PerforationPressureGradient(backend);
            props.WellIndex = WellIndex(backend);
            props.FacilityWellMapping = FacilityWellMapping(backend);
        end
    end
end