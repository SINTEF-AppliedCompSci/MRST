classdef FluxDiscretization < PropertyFunctions
    properties
        PermeabilityPotentialGradient % K * (grad(p) + rho g dz)
        PressureGradient % Gradient of phase pressures
        GravityPotentialDifference % rho * g * dz term
        PhasePotentialDifference % (grad p_alpha + dpdz)
        PhaseFlux % Phase volumetric fluxes
        FaceMobility % Phase mobility on face
        FaceComponentMobility % Composition * mobility on face
        PhaseUpwindFlag
        ComponentTotalFlux % Total mass flux for each component
        ComponentPhaseFlux % Phase fluxes for each component
        Transmissibility
    end

    properties (Access = protected)
        FlowStateBuilder
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
            
            % Face values - typically upwinded
            props.FaceComponentMobility = FaceComponentMobility(backend, upstr);
            props.FaceMobility = FaceMobility(backend, upstr);
            % 
            props.ComponentPhaseFlux = ComponentPhaseFlux(backend);
            props.ComponentTotalFlux = ComponentTotalFlux(backend);
            props.PhaseFlux = PhaseFlux(backend);

            % Flow discretizer
            props.FlowStateBuilder = ImplicitFlowStateBuilder();
            
            % Define storage
            props.structName = 'FluxProps';
        end

        function fd = setFlowStateBuilder(fd, fb)
            fd.FlowStateBuilder = fb;
        end
        
        function fb = getFlowStateBuilder(fd)
            fb = fd.FlowStateBuilder;
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
            flowState = fd.buildFlowState(model, state, state0, dt);
            v = model.getProps(flowState, 'ComponentTotalFlux');
            for c = 1:ncomp
                acc{c} = (mass{c} - mass0{c})./dt;
            end
        end

        function flowState = buildFlowState(fd, model, state, state0, dt)
            flowState = fd.FlowStateBuilder.build(fd, model, state, state0, dt);
        end
        
        function dt = getMaximumTimestep(fd, model, state, state0, dt, forces)
            dt = fd.FlowStateBuilder.getMaximumTimestep(fd, model, state, state0, dt, forces);
        end
        
        function [fd, state] = prepareTimestep(fd, model, state, state0, dt, drivingForces)
            % Called before each time-step solve
            [fd.FlowStateBuilder, state] = fd.FlowStateBuilder.prepareTimestep(fd, model, state, state0, dt, drivingForces);
        end
    end
end