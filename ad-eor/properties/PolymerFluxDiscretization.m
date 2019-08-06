classdef PolymerFluxDiscretization < FluxDiscretization
    properties
        PolymerPhaseFlux % Polymer phase volumetric fluxes
        FaceConcentration % Polymer or surfactant concentration on face
    end

    methods
        function props = PolymerFluxDiscretization(model)
            props = props@FluxDiscretization(model);
            upstr = UpwindFunctionWrapperDiscretization(model);
            props.PolymerPhaseFlux = PolymerPhaseFlux(model);            
            props.FaceConcentration = FaceConcentration(model, upstr);
        end
    end

end