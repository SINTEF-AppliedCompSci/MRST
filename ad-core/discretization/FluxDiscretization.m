classdef FluxDiscretization < PropertyFunctions
    properties
        PermeabilityPotentialGradient % K * (grad(p) + rho g dz)
        PressureGradient % Gradient of phase pressures
        GravityPotentialDifference % rho * g * dz term
        PhasePotentialDifference % (grad p_alpha + dpdz)
        PhaseFlux % Phase volumetric fluxes
        PhaseUpwindFlag
        ComponentFlux
        ComponentMass
        Transmissibility
    end
    
    methods
        function props = FluxDiscretization(model)
            backend = model.AutoDiffBackend;
            upstr = UpstreamFunctionWrapper(model.operators.faceUpstr);
            tpfa = TwoPointFluxApproximation(model);

            props@PropertyFunctions();
            % Darcy flux
            props.Transmissibility = Transmissibility(backend);
            props.PermeabilityPotentialGradient = PermeabilityPotentialGradient(backend, tpfa);
            props.PressureGradient = PressureGradient(backend);
            props.GravityPotentialDifference = GravityPotentialDifference(backend);
            
            % Phase flux
            props.PhasePotentialDifference = PhasePotentialDifference(backend);
            props.PhaseUpwindFlag = PhaseUpwindFlag(backend);
            % Components
            props.ComponentMass = ComponentMass(backend);
            
            % Fluxes - these are upwinded properties
            props.ComponentFlux = ComponentFlux(backend, upstr);
            props.PhaseFlux = PhaseFlux(backend, upstr);

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