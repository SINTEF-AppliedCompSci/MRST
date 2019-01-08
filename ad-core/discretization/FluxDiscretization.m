classdef FluxDiscretization < PropertyFunctions
    properties
        PermeabilityPotentialGradient % K * (grad(p) + rho g dz)
        PressureGradient % Gradient of phase pressures
        GravityPotentialDifference % rho * g * dz term
        PhasePotentialDifference % (grad p_alpha + dpdz)
        FaceMobility % Mobility on the face (for advective transport)
        PhaseFlux % Phase volumetric fluxes
        ComponentMassFlux
    end
    
    methods
        function props = FluxDiscretization(model)
            backend = model.AutoDiffBackend;
            upstr = UpstreamFunctionWrapper(model.operators.faceUpstr);
            tpfa = TwoPointFluxApproximation(model);

            props@PropertyFunctions();
            props.PermeabilityPotentialGradient = PermeabilityPotentialGradient(backend, tpfa);
            props.PressureGradient = PressureGradient(backend);
            props.GravityPotentialDifference = GravityPotentialDifference(backend);
            props.PhasePotentialDifference = PhasePotentialDifference(backend);
            props.FaceMobility = FaceMobility(backend, upstr);
            props.PhaseFlux = PhaseFlux(backend);
            props.ComponentMassFlux = [];

            
            % Define storage
            props.structName = 'FluxProps';
        end

        function state = evaluateProperty(props, model, state, name)
            switch name

            end
            state = evaluateProperty@PropertyFunctions(props, model, state, name);
        end
        
        function v = faceUpstreamPhase(fd, state, flag, cellvalue)
            v = fd.PhaseUpstreamDiscretization.faceUpstream(state, flag, cellvalue);
        end
        
        function v = faceUpstreamComponent(fd, state, flag, cellvalue)
            v = fd.ComponentUpstreamDiscretization.faceUpstream(state, flag, cellvalue);
        end
    end
end