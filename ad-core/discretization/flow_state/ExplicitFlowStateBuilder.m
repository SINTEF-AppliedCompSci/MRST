classdef ExplicitFlowStateBuilder < FlowStateBuilder
    properties
        saturationCFL = 1;
        compositionCFL = 1;
        explicitFlowProps = {'FaceMobility', 'FaceComponentMobility',...
                             'CapillaryPressure', 'GravityPotentialDifference'};
        implicitFlowProps = {'PressureGradient'};
        initialStep = 1*day;
    end
    
    methods
        function dt_max = getMaximumTimestep(fsb, fd, model, state, state0, dt, forces)
            if ~isfield(state, 'flux')
                dt_max = fsb.initialStep;
                return;
            end
            
            if isa(model, 'ThreePhaseCompositionalModel') && isfinite(fsb.compositionCFL)
                
            end
            cfl_s = estimateSaturationCFL(model, state, 1, 'forces', forces);
            dt_max = 1./max(cfl_s);
        end
        
        function flowState = build(builder, fd, model, state, state0, dt)
            % Hybridize state
            % Get implicit props to ensure they are cached
            model.getProps(state, builder.implicitFlowProps{:});
            flowState = state;
            props = builder.explicitFlowProps;
            name = fd.getPropertyContainerName();
            for i = 1:numel(props)
                if isfield(state0, name)
                    % Remove cached entries
                    if ~isempty(state0.(name).(prop))
                        state0.(name).(prop) = [];
                    end
                end
                prop = props{i};
                f = model.getProps(state0, prop);
                flowState.(name).(prop) = f;
            end
        end
        
        function [builder, state] = prepareTimestep(builder, fd, model, state, state0, dt, drivingForces)
            assert(model.outputFluxes, 'model.outputFluxes must be true for explicit solver');
        end
    end
end
