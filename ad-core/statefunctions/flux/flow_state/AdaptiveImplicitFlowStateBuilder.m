classdef AdaptiveImplicitFlowStateBuilder < ExplicitFlowStateBuilder
    % AIM - adaptive implicit flow state builder. Takes certain cells to be
    % implicit based on estimated CFl.
    properties
        firstStepImplicit = true; % The first step is taken implicit if fluxes are missing
        implicitWells = true; % Well cells are always implicit
    end
    
    methods
        function builder = AdaptiveImplicitFlowStateBuilder(varargin)
            builder@ExplicitFlowStateBuilder(varargin{:});
            builder.explicitFlux = {}; % Not supported
            builder.explicitFlow = {'Mobility', 'ComponentMobility'};
        end
        
        function dt_max = getMaximumTimestep(fsb, fd, model, state, state0, dt, forces)
            % No time-step limit for AIM
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
            [groups, names] = builder.getExplicitGroups(model);
            for groupNo = 1:numel(groups)
                grp = groups{groupNo};
                props = names{groupNo};
                name = grp.getStateFunctionContainerName();
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
                props = builder.explicitProps;
                for i = 1:numel(props)
                    p = props{i};
                    v0 = model.getProp(state0, p);
                    v = model.getProp(state, p);
                    if iscell(v)
                        for j = 1:numel(v)
                            v{j} = v{j}.*implicit;
                            if iscell(v0)
                                vj0 = v0{j};
                            else
                                vj0 = v0(:, j);
                            end
                            v{j} = v{j} + explicit.*vj0;
                        end
                    else
                        v = v.*implicit + explicit.*v0;
                    end
                    flowState = model.setProp(flowState, p, v);
                end
            end
        end
        
        function [builder, state] = prepareTimestep(builder, fd, model, state, state0, dt, drivingForces)
            [builder, state] = prepareTimestep@ExplicitFlowStateBuilder(builder, fd, model, state, state0, dt, drivingForces);
            isFirst = builder.isFirstTimeStep(state);
            if isFirst
                fs = builder.firstStepImplicit;
                [impl_sat, impl_comp] = deal(repmat(fs, model.G.cells.num, 1));
            else
                iflow = builder.useInflowForEstimate;
                cfl = estimateSaturationCFL(model, state, dt, 'forces', drivingForces, 'useInflow', iflow);
                impl_sat = cfl >= builder.saturationCFL;
                
                cfl_c = estimateCompositionCFL(model, state, dt, 'forces', drivingForces, 'useInflow', iflow);
                impl_comp = max(cfl_c, [], 2) >= builder.compositionCFL;
            end
            implicit = impl_sat | impl_comp;
            if builder.implicitWells && ~isempty(drivingForces.W)
                % Well cells are taken implicitly
                wc = vertcat(drivingForces.W.cells);
                implicit(wc) = true;
            else
                wc = [];
            end
            nc = numel(implicit);
            ni_s = sum(impl_sat);
            ni_c = sum(impl_comp);
            ni_w = numel(wc);
            ni = sum(implicit);
            dispif(builder.verbose, ...
                ['Adaptive implicit: %d of %d cells are implicit (%2.2f%%).\n', ...
                '%d limited by composition, %d limited by saturation, %d belong to wells.\n'], ...
                ni, nc, 100*ni/nc, ni_c, ni_s, ni_w);
            state.implicit = implicit;
            assert(isempty(builder.explicitFlux), ...
                'AIM does not support flux terms. You should specify PVT or Flow propertes instead.');
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
