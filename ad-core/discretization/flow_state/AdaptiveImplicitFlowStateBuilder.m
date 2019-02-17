classdef AdaptiveImplicitFlowStateBuilder < ExplicitFlowStateBuilder
    properties
        
    end
    
    methods
        function builder = AdaptiveImplicitFlowStateBuilder(varargin)
            builder@ExplicitFlowStateBuilder(varargin{:});
            builder.explicitFlowProps = {'Mobility', 'ComponentMobility'};
            builder.implicitFlowProps = {};
        end
        
        function dt_max = getMaximumTimestep(fsb, fd, model, state, state0, dt, forces)
            dt_max = inf;
        end
        
        function flowState = build(builder, fd, model, state, state0, dt)
            % Hybridize state
            flowState = state;
            implicit = state.implicit;
            if all(implicit)
                return
            end
            explicit = ~implicit;
            props = builder.explicitFlowProps;
            
            fp = model.FlowPropertyFunctions;
            name = fp.getPropertyContainerName();
            for i = 1:numel(props)
                prop = props{i};
                if isfield(state0, name)
                    % Remove cached entries
                    if ~isempty(state0.(name).(prop))
                        state0.(name).(prop) = [];
                    end
                end
                X = model.getProps(state, prop);
                X0 = model.getProps(state0, prop);
                
                X_hyb = X;
                if iscell(X_hyb)
                    for j = 1:numel(X_hyb)
                        if ~isempty(X_hyb{j})
                            X_hyb{j} = implicit.*X{j} + explicit.*X0{j};
                        end
                    end
                else
                    X_hyb = implicit.*X + explicit.*X0;
                end
                flowState.(name).(prop) = X_hyb;
            end
        end
        
        function [builder, state] = prepareTimestep(builder, fd, model, state, state0, dt, drivingForces)
            [builder, state] = prepareTimestep@ExplicitFlowStateBuilder(builder, fd, model, state, state0, dt, drivingForces);
            if isfield(state, 'flux')
                cfl = estimateSaturationCFL(model, state, dt, 'forces', drivingForces);
                implicit = cfl >= builder.saturationCFL;
            else
                implicit = true(model.G.cells.num, 1);
            end
            nc = numel(implicit);
            ni = sum(implicit);
            dispif(builder.verbose, 'Adaptive implicit: %d of %d cells are implicit (%2.2f%%)\n', ni, nc, 100*ni/nc);
            state.implicit = implicit;
        end
    end
end
