classdef FluxDiscretizationDG < FluxDiscretization
   
    properties
        parent
        ComponentTotalVelocity
        ComponentPhaseVelocity
        TotalVelocity
    end
    
    methods
        function fd = FluxDiscretizationDG(model)
 
            fd = fd@FluxDiscretization(model);
            
            fd = fd.setStateFunction('PhaseFlux', PhaseFluxFixedTotalVelocity(model));
            fd = fd.setStateFunction('PhaseUpwindFlag', PhasePotentialUpwindFlag(model));
            fd = fd.setStateFunction('ComponentPhaseFlux', ComponentPhaseFluxFractionalFlow(model));
            
            fd = fd.setStateFunction('PhaseInterfacePressureDifferences', PhaseInterfacePressureDifferences(model));
            fd = fd.setStateFunction('TotalFlux', FixedTotalFluxDG(model));
            fd = fd.setStateFunction('FaceTotalMobility', FaceTotalMobility(model));
            
            fd = fd.setStateFunction('Transmissibility', TransmissibilityDG(model));
            fd = fd.setStateFunction('GravityPotentialDifference', GravityPotentialDifferenceDG(model));
            fd = fd.setStateFunction('PhaseUpwindFlag', PhasePotentialUpwindFlagDG(model));
            
            fd.ComponentTotalVelocity = ComponentTotalVelocityDG(model);
            fd.ComponentPhaseVelocity = ComponentPhaseVelocityFractionalFlowDG(model);
            fd.TotalVelocity = FixedTotalVelocityDG(model);
            
            fd = fd.setFlowStateBuilder(FlowStateBuilderDG);
        end
        
        function [acc, v, vc, names, types] = componentConservationEquations(fd, model, state, state0, dt)
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
            mass = model.getProps(state.cellStateDG, 'ComponentTotalMass');
            mass0 = model.getProps(state0.cellStateDG, 'ComponentTotalMass');
            flowState = fd.buildFlowState(model, state, state0, dt);
            v  = model.getProps(flowState, 'ComponentTotalFlux');
            flowState.FlowProps = state.cellStateDG.FlowProps;
            flowState.cells = state.cellStateDG.cells;
            vc = model.getProps(flowState, 'ComponentTotalVelocity');
            for c = 1:ncomp
                acc{c} = (mass{c} - mass0{c})./dt;
            end
        end
        
    end
    
end