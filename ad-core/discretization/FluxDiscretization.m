classdef FluxDiscretization < PropertyFunctions
    properties
        PermeabilityPotentialGradient % K * (grad(p) + rho g dz)
        PressureGradient % Gradient of phase pressures
        GravityPotentialDifference % rho * g * dz term
        PhasePotentialDifference % (grad p_alpha + dpdz)
        PhaseFlux % Phase volumetric fluxes
        PhaseUpwindFlag
        ComponentTotalFlux % Total mass flux for each component
        ComponentPhaseFlux % Phase fluxes for each component
        Transmissibility
    end

    properties (Access = private)
        
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
            
            % Fluxes - these are upwinded properties
            props.ComponentPhaseFlux = ComponentPhaseFlux(backend, upstr);
            props.ComponentTotalFlux = ComponentTotalFlux(backend);
            props.PhaseFlux = PhaseFlux(backend, upstr);

            % Define storage
            props.structName = 'FluxProps';
        end

        function state = evaluateProperty(props, model, state, name)
            switch name

            end
            state = evaluateProperty@PropertyFunctions(props, model, state, name);
        end
        
        function [acc, v, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            ncomp = model.getNumberOfComponents;
            [acc, types] = deal(cell(1, ncomp));
            names = model.getComponentNames();
            [types{:}] = deal('cell');
            mass = model.getProps(state, 'ComponentTotalMass');
            mass0 = model.getProps(state0, 'ComponentTotalMass');
            v = model.getProps(state, 'ComponentTotalFlux');
            for c = 1:ncomp
                acc{c} = (mass{c} - mass0{c})./dt;
            end
        end
    end
end