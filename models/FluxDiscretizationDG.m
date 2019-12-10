classdef FluxDiscretizationDG < FluxDiscretization
   
    properties
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
            % Get mass terms
            model = fd.expandRegions(model, 'cells');
            mass  = model.getProps(state.cellStateDG, 'ComponentTotalMass');
            mass0 = model.getProps(state0.cellStateDG, 'ComponentTotalMass');
            % Get boundary flux term
            flowState = fd.buildFlowState(model, state, state0, dt);
            model     = fd.expandRegions(model, 'faces');
            v         = model.getProps(flowState, 'ComponentTotalFlux');
            % Get cell flux term
            state.cellStateDG.FluxProps = state.faceStateDG.FluxProps;
            model = fd.expandRegions(model, 'cells');
            vc = model.getProps(state.cellStateDG, 'ComponentTotalVelocity');
            for c = 1:ncomp
                acc{c} = (mass{c} - mass0{c})./dt;
            end
        end
        
        function model = expandRegions(fd, model, type, elements)
            
            if nargin < 4
                elements = [];
            end
            switch type
                case 'cells'
                    cells = elements;
                    if isempty(cells)
                        [~, ~, cells] = model.discretization.getCubature((1:model.G.cells.num)', 'cell');
                    end
                case 'faces'
                    faces = elements;
                    if isempty(faces)
                        [~, ~, ~, faces] = model.discretization.getCubature(find(model.operators.internalConn), 'face');
                    end
                    cells = [model.G.faces.neighbors(faces,1); model.G.faces.neighbors(faces,2)];
            end
            
            fpprops = model.FlowPropertyFunctions;
            fpprops = expandPropsRegions(fpprops, cells);
            model.FlowPropertyFunctions = fpprops;
            
            fdprops = model.FluxDiscretization;
            fdprops = expandPropsRegions(fdprops, cells);
            model.FluxDiscretization = fdprops;
            
        end
        
    end
    
end

function props = expandPropsRegions(props, cells)
    names = props.getNamesOfStateFunctions();
    for j = 1:numel(names)
        if isprop(props, names{j})
            p = props.(names{j});
            if isprop(p, 'regions') && ~isempty(p.regions)
                p.regions = p.regions(cells);
                props.(names{j}) = p;
            end
        end
    end
end