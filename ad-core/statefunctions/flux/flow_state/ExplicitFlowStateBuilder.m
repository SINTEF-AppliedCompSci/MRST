classdef ExplicitFlowStateBuilder < FlowStateBuilder
    properties
        saturationCFL = 0.9;
        compositionCFL = 0.9;
        explicitFluxProps = {'FaceMobility', 'FaceComponentMobility',...
                             'GravityPotentialDifference'};
        implicitFluxProps = {'PressureGradient'};
        initialStep = 1*day;
    end
    
    methods
        function dt_max = getMaximumTimestep(fsb, fd, model, state, state0, dt, forces)
            if ~isfield(state, 'flux')
                dt_max = fsb.initialStep;
                return;
            end
            cfl_s = estimateSaturationCFL(model, state, 1/fsb.saturationCFL, 'forces', forces);
            cfl_c = estimateCompositionCFL(model, state, 1/fsb.compositionCFL, 'forces', forces);

            dt_max = 1./max([cfl_s, cfl_c], [], 2);
            dt_max = min(dt_max);
        end
        
        function flowState = build(builder, fd, model, state, state0, dt)
            % Hybridize state
            % Get implicit props to ensure they are cached
            model.getProps(state, builder.implicitFluxProps{:});
            flowState = state;
            props = builder.explicitFluxProps;
            name = fd.getStateFunctionContainerName();
            if ~isfield(state0, name)
                % Ensure that property containers exist
                state0 = model.initStateFunctionContainers(state0);
                state0 = value(state0);
            end
            if ~isfield(flowState, name)
                flowState = model.initStateFunctionContainers(flowState);
            end
            for i = 1:numel(props)
                prop = props{i};
                % Remove cached entries
                if ~isempty(state0.(name).(prop))
                    state0.(name).(prop) = [];
                end
                f = model.getProps(state0, prop);
                flowState.(name).(prop) = f;
            end
        end
        
        function [builder, state] = prepareTimestep(builder, fd, model, state, state0, dt, drivingForces)
            assert(model.outputFluxes, 'model.outputFluxes must be true for explicit solver');
        end
    end
end
