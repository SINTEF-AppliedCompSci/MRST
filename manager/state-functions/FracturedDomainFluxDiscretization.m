classdef FracturedDomainFluxDiscretization < StateFunctionGrouping
    % Function grouping for Darcy-type flux discretization. The defaults
    % gives a industry-standard single-point upwind scheme with a two-point
    % flux discretization which emphasizes robustness and efficiency.
    properties
        PermeabilityPotentialGradient % K * (grad(p) + rho g dz)
        PressureGradient % Gradient of phase pressures
        GravityPotentialDifference % rho * g * dz term
        PhasePotentialDifference % (grad p_alpha + dpdz)
        PhaseFlux % Phase volumetric fluxes
        FaceMobility % Phase mobility on face
        FaceComponentMobility % Composition * mobility on face
        PhaseUpwindFlag % Upwind flag for each phase
        ComponentTotalFlux % Total mass flux for each component
        ComponentPhaseFlux % Phase fluxes for each component
        Transmissibility % Face-based transmissibility
    end

    properties (Access = protected)
        FlowStateBuilder
    end
    
    methods
        function props = FracturedDomainFluxDiscretization(model)
            upstr = UpwindFunctionWrapperDiscretization(model);
            tpfa = TwoPointFluxApproximation(model);

            props@StateFunctionGrouping();
            % Darcy flux
            if ~isempty(model.inputdata) && isfield(model.inputdata.SOLUTION, 'THPRES')
                trans = ThresholdedTransmissibility(model, model);
            else
                trans = Transmissibility(model);
            end
            props.Transmissibility = trans;
            props.PermeabilityPotentialGradient = PermeabilityPotentialGradient(model, tpfa);
            props.PressureGradient = PressureGradient(model);
            props.GravityPotentialDifference = FracturedDomainGravityPotentialDifference(model);
            
            % Phase flux
            props.PhasePotentialDifference = PhasePotentialDifference(model);
            props.PhaseUpwindFlag = PhaseUpwindFlag(model);
            
            % Face values - typically upwinded
            props.FaceComponentMobility = FaceComponentMobility(model, upstr);
            props.FaceMobility = FracturedDomainFaceMobility(model);
            % 
            props.ComponentPhaseFlux = ComponentPhaseFlux(model);
            props.ComponentTotalFlux = ComponentTotalFlux(model);
            props.PhaseFlux = FracturedDomainPhaseFlux(model);

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
        
        function [acc, v, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            % Compute discretized conservation equations in the interior of the domain.
            % REQUIRED PARAMETERS:
            %   model  - ReservoirModel derived class
            %   state  - State corresponding to the given model at the end
            %            of the current time-integration interval.
            %   state0 - State at previous time-step.
            %   dt     - Time-step length.
            % RETURNS:
            %   acc    - Cell array of component accumulation terms
            %            discretized with a first order finite-difference
            %            scheme.
            %   v      - Cell array of component fluxes in the interior of
            %            the domain.
            %   names  - The names of each equation (corresponds to
            %            component names.
            %   types  - Types of each equation (defaults to 'cell', to
            %            indicate cell-wise equations on the whole domain)
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