classdef ExplicitFlowStateBuilder < FlowStateBuilder
    properties
        saturationCFL = 0.9; % Target saturation CFL. Should be <= 1 for stability.
        compositionCFL = 0.9; % Target composition CFL. Should be <= 1 for stability.
        explicitFlux = {'FaceMobility', 'FaceComponentMobility',...
                        'GravityPotentialDifference'}; % StateFunctions that should always be explicit
        explicitFlow = {}; % Flow properties
        explicitPVT = {}; % PVT properties
        explicitProps = {}; % setProp capable properties that should be explicit
        initialStep = 1*day;% Timestep used if no fluxes are present
        useInflowForEstimate = false;
        changeTimestep = true; % Estimate stable step. Can be disabled if you know your steps are small enough.
    end
    
    methods
        function fsb = ExplicitFlowStateBuilder(varargin)
            fsb@FlowStateBuilder(varargin{:});
        end
        
        function dt_max = getMaximumTimestep(fsb, fd, model, state, state0, dt, forces)
            if ~fsb.changeTimestep
                dt_max = inf;
                return
            end
            if fsb.isFirstTimeStep(state)
                dt_max = fsb.initialStep;
                return;
            end
            % Remove any cached properties
            state = model.reduceState(state, true);
            iflow = fsb.useInflowForEstimate;
            cfl_s = estimateSaturationCFL(model, state, 1/fsb.saturationCFL, 'forces', forces, 'useInflow', iflow);
            cfl_c = estimateCompositionCFL(model, state, 1/fsb.compositionCFL, 'forces', forces, 'useInflow', iflow);
            dt_max_z = min(1./max(cfl_c, [], 2));
            dt_max_s = min(1./max(cfl_s, [], 2));
            dt_max = min(dt_max_z, dt_max_s);
            if fsb.verbose
                if dt_max < dt
                    if dt_max_z > dt_max_s
                        s = 'composition';
                    else
                        s = 'saturation';
                    end
                    fprintf('Time-step limited by %s CFL: %s reduced to %s (%1.2f%% reduction)\n', ...
                        s, formatTimeRange(dt, 0), formatTimeRange(dt_max, 0), 100*(dt - dt_max)/dt);
                end
            end
        end
        
        function flowState = build(builder, fd, model, state, state0, dt)
            % Hybridize state
            % The base state is the implicit. Other functions are then
            % assigned.
            flowState = state;
            % First insert any properties that can be set with setProp
            props = builder.explicitProps;
            for i = 1:numel(props)
                p = props{i};
                flowState = model.setProp(flowState, p, model.getProp(state0, p));
            end
            % Then we make sure that the explicit flux state functions are
            % set from the explicit state.
            [groups, names] = builder.getExplicitGroups(model);
            for groupNo = 1:numel(groups)
                grp = groups{groupNo};
                props = names{groupNo};
                name = grp.getStateFunctionContainerName();
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
        end
        
        function [groups, names] = getExplicitGroups(builder, model)
            groups = {};
            names = {};
            pvt = builder.explicitPVT;
            if ~isempty(pvt)
                groups{end+1} = model.PVTPropertyFunctions;
                names{end+1} = pvt;
            end
            flow = builder.explicitFlow;
            if ~isempty(flow)
                groups{end+1} = model.FlowPropertyFunctions;
                names{end+1} = flow;
            end
            flux = builder.explicitFlux;
            if ~isempty(flux)
                groups{end+1} = model.FlowDiscretization;
                names{end+1} = flux;
            end
        end
        
        function isFirst = isFirstTimeStep(builder, state)
            isFirst = ~isfield(state, 'flux') || ~any(max(abs(state.flux)));
        end
        
        function [builder, state] = prepareTimestep(builder, fd, model, state, state0, dt, drivingForces)
            assert(model.outputFluxes, 'model.outputFluxes must be true for explicit solver');
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
