classdef FluxDiscretization < PropertyFunctions
    properties
        PermPotentialGradient
        PotentialGradient
        FaceGravityPotential
        FaceMobility
        PhaseFlux
        ComponentMassFlux
    end
    
    methods
        function props = FluxDiscretization(model)
            props@PropertyFunctions();
            % Saturation properties
            props.FaceGravityPotential = FaceGravityPotential(model.AutoDiffBackend);
            props.PotentialGradient = PotentialGradient(model.AutoDiffBackend);
            % Define storage
            props.structName = 'FluxProps';
        end

        function state = evaluateProperty(props, model, state, name)
            switch name

            end
            state = evaluateProperty@PropertyFunctions(props, model, state, name);
        end
    end
end