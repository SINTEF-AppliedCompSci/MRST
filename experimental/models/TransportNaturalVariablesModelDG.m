classdef TransportNaturalVariablesModelDG < TransportNaturalVariablesModel
    properties
        disc
        tryMaxDegree
    end
    
    methods
        function model = TransportNaturalVariablesModelDG(G, rock, fluid, compFluid, varargin)
            
            model = model@TransportNaturalVariablesModel(G, rock, fluid, compFluid, varargin{:});
            model.tryMaxDegree = true;
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = transportEquationsNaturalVarsDG(state0, state, model, dt, ...
                            drivingForces, varargin{:});
        end
        
        % ----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name)
            % Map variables to state field.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.getVariableField`
            switch(lower(name))
                case {'swdof'}
                    index = model.satVarIndex('sWdof');
                    fn = 'sdof';
                case {'sodof'}
                    index = model.satVarIndex('sOdof');
                    fn = 'sdof';
                case {'sgdof'}
                    index = model.satVarIndex('sGdof');
                    fn = 'sdof';
                case {'sdof'}
                    index = ':';
                    fn = 'sdof';
                case {'xdof'}
                    index = ':';
                    fn = 'xdof';
                case {'ydof'}
                    index = ':';
                    fn = 'ydof';
                case {'componentsdof'}
                    index = ':';
                    fn = 'componentsdof';
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@TransportNaturalVariablesModel(model, name);
            end
        end
        
        % ----------------------------------------------------------------%
        function vars = getSaturationVarNames(model)
            vars = {'sWdof', 'sOdof', 'sGdof'};
            ph = model.getActivePhases();
            vars = vars(ph);
        end
        
        % ----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            
            state0 = state;
            dx0    = dx;
            vars   = problem.primaryVariables;
            if model.allowLargeSaturations
                dsMax = max(sum(state.sdof, 2), 1).*model.dsMaxAbs;
            else
                dsMax = model.dsMaxAbs;
            end
            [dxMax, dyMax] = deal(model.dzMaxAbs);
            
             % Get saturation and composition updates
            [ds, deltax, deltay, dL, dx, vars] = model.getHyperbolicUpdates(problem, dx, vars, state);

            w_sat = min(dsMax./max(abs(ds), [], 2), 1);
            w_x = min(dxMax./max(abs(deltax), [], 2), 1);
            w_y = min(dyMax./max(abs(deltay), [], 2), 1);
            w_xy = min(w_x, w_y);
                        
            % Minimum relax factor
            w = min(w_sat, w_xy);
            % Update sat
            capunit = @(x) min(max(x, 0), 1);
            if any(ds(:) ~= 0)
                if model.water
                    ds_relax = ds;
                    ds_relax(:, 1) = ds(:, 1).*min(dsMax./max(abs(ds(:, 1)), [], 2), 1);
                    ds_relax(:, 2:end) = ds(:, 2:end).*w;
                else
                    ds_relax = bsxfun(@times, w, ds);
                end
%                 sdof_uncap = state.sdof + ds_relax;
                state.sdof = state.sdof + ds_relax;
                state.s = state.s.*0;
                for phNo = 1:size(state.s,2)
                    state.s(:,phNo) = model.disc.getCellMean(state.sdof(:,phNo), state);
                end
                s_uncap = state.s;
                if model.allowLargeSaturations
                    bads = any(state.s < 0 - model.disc.meanTolerance, 2);
%                     if any(bads)
%                         state.s(bads,:) = max(state.s(bads,:), 0);
%                         state = dgLimiter(model.disc, state, bad, 's', 'kill');
%                     end
                else
                    bads = any((state.s > 1 + model.disc.meanTolerance) ...
                             | (state.s < 0 - model.disc.meanTolerance), 2);
                
%                     if any(bad)
%                         state.s(bad, :) = min(state.s(bad, :), 1);
%                         state.s(bad, :) = max(state.s(bad, :), 0);
%                         state.s(bad, :) = bsxfun(@rdivide, state.s(bad, :), ...
%                                                        sum(state.s(bad, :), 2));
% 
%                         state = dgLimiter(model.disc, state, bad, 's', 'kill');
%                     end
                end
            else
                s_uncap = state.s;
            end
            
            % Update dx, dy
            if any(deltax(:) ~= 0) || any(deltay(:) ~= 0)
                state.xdof = max(state.xdof + bsxfun(@times, deltax, w), 0);
                state.ydof = max(state.ydof + bsxfun(@times, deltay, w), 0);
                
                for cNo = 1:size(state.xdof,2)
                    state.x(:,cNo) = model.disc.getCellMean(state.xdof(:,cNo), state);
                    state.y(:,cNo) = model.disc.getCellMean(state.ydof(:,cNo), state);
                end
                
                badx = any((state.x > 1 + model.disc.meanTolerance) ...
                         | (state.x < 0 - model.disc.meanTolerance), 2);
                bady = any((state.y > 1 + model.disc.meanTolerance) ...
                         | (state.y < 0 - model.disc.meanTolerance), 2);
                
%                 state.xdof = capunit(state.xdof);
%                 state.xdof = bsxfun(@rdivide, state.xdof, sum(state.xdof, 2));
% 
%                 state.y = capunit(state.y);
%                 state.y = bsxfun(@rdivide, state.y, sum(state.y, 2));
                
                xyUpdated = true;
            else
                xyUpdated = false;
            end
            
            bad = bads | badx | bady;
            if any(bad)
                state.s(bad,:) = max(state.s(bad,:), 0);
                degree = state.degree;
                state = dgLimiter(model.disc, state, bad, 's', 'kill');
                state.degree = degree;
                state = model.disc.updateDofPos(state);
                state.x(bad,:) = capunit(state.x(bad,:));
                state.x = bsxfun(@rdivide, state.x, sum(state.x, 2));
                state = dgLimiter(model.disc, state, bad, 'x', 'kill', 'plot', model.disc.plotLimiterProgress);
                state.degree = degree;
                state = model.disc.updateDofPos(state);
                state.y(bad,:) = capunit(state.y(bad,:));
                state.y = bsxfun(@rdivide, state.y, sum(state.y, 2));
                state = dgLimiter(model.disc, state, bad, 'y', 'kill');
%                 state = model.disc.updateDofPos(state);
            end
            
            if ~isempty(dL)
                state.L = capunit(state.L + w.*dL);
            end
            
            problem.primaryVariables = vars;
            [state, report] = updateState@ThreePhaseCompositionalModel(model, state, problem, dx, drivingForces);
            statePrev = state;
            state  = model.flashPhases(state, state0, s_uncap, xyUpdated, problem.iterationNo);
            
            ds = state.s./statePrev.s;
            ds(isnan(ds)) = 1;
            state.sdof = state.sdof.*rldecode(ds, state.nDof, 1);
            for phNo = 1:size(state.s,2)
                ix = model.disc.getDofIx(state, 1, isinf(ds(:,phNo)));
                state.sdof(ix,phNo) = state.s(isinf(ds(:,phNo)), phNo);
            end
            for phNo = 1:size(state.sdof,2)
                state.s(:,phNo) = model.disc.getCellMean(state.sdof(:,phNo), state);
            end
            
            dx = state.x./statePrev.x;
            dx(isnan(dx)) = 1;
            state.xdof = state.xdof.*rldecode(dx, state.nDof, 1);
            for cNo = 1:size(state.x,2)
                ix = model.disc.getDofIx(state, 1, isinf(dx(:,cNo)));
                state.xdof(ix,cNo) = state.x(isinf(dx(:,cNo)), cNo);
            end
            
            dy = state.y./statePrev.y;
            dy(isnan(dy)) = 1;
            state.ydof = state.ydof.*rldecode(dy, state.nDof, 1);
            for cNo = 1:size(state.y,2)
                ix = model.disc.getDofIx(state, 1, isinf(dy(:,cNo)));
                state.ydof(ix,cNo) = state.y(isinf(dy(:,cNo)), cNo);
            end
            
            for cNo = 1:size(state.xdof,2)
                state.x(:,cNo) = model.disc.getCellMean(state.xdof(:,cNo), state);
                state.y(:,cNo) = model.disc.getCellMean(state.ydof(:,cNo), state);
            end
            
            zdof = bsxfun(@times, state.xdof, rldecode(state.L, state.nDof, 1)) ...
                 + bsxfun(@times, state.ydof,  rldecode(1-state.L, state.nDof, 1));
            state.componentsdof = zdof;
            
            if model.water
                sT = sum(state.s, 2);
                space = (sT - state.s(:, 1))./sT;
            else
                space = 1;
            end
            
            state.dz = max(max(abs(bsxfun(@times, deltax, space)), bsxfun(@times, deltay, space)));
            dz2 = max(abs(bsxfun(@times, state.components - state0.components, space)));
            state.dz = max(state.dz, dz2);
            
            
        end

        
        % ----------------------------------------------------------------%
        function [ds, deltax, deltay, dL, dx, vars] = getHyperbolicUpdates(model, problem, dx, vars, state)
            [ds, dx, vars] = model.getSaturationIncrements(dx, vars, state);
            [deltax, deltay, dx, vars] = model.getPhaseCompositionIncrements(dx, vars, state);
            
%             [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
            pureLiquid = state.pureLiquid;
            mixvar = strcmpi(vars, 'sGsO');
            if any(mixvar)
                % Sequential implicit stuff
                dMix = dx{mixvar};
                ix = model.disc.getDofIx(state, Inf, ~pureLiquid);
                ds(ix, end) = dMix(ix);
                ix = model.disc.getDofIx(state, Inf, pureLiquid);
                ds(ix, end-1) = dMix(ix);
                [vars, removed] = model.stripVars(vars, 'sGsO');
                dx = dx(~removed);
            end
            dL = [];
        end
        
        % ----------------------------------------------------------------%
        function [ds, dx, vars] = getSaturationIncrements(model, dx, vars, state)
%             [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
            pureLiquid = state.pureLiquid;
            pureVapor  = state.pureVapor;
            twoPhase   = state.twoPhase;
            [dsG, dx, vars] = getSatUpdateInternal(model, 'satg', dx, vars, twoPhase, state);
            [dsO, dx, vars] = getSatUpdateInternal(model, 'sato', dx, vars, twoPhase, state);
            [dsW, dx, vars] = getSatUpdateInternal(model, 'satw', dx, vars, twoPhase, state);
            if model.water
                if ~any(strcmpi(vars, 'sGsO'))
                    dsO(pureLiquid) = -dsW(pureLiquid);
                    dsG(pureVapor) = -dsW(pureVapor);
                end
                ds = [dsW, dsO, dsG];
            else
                ds = [dsO, dsG];
            end
        end

        % ----------------------------------------------------------------%
        function [dx, dy, increments, vars] = getPhaseCompositionIncrements(model, increments, vars, state)
%             [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
            pureLiquid = state.pureLiquid;
            pureVapor  = state.pureVapor;
            twoPhase   = state.twoPhase;
            pureVaporIx = model.disc.getDofIx(state, Inf, pureVapor);
            twoPhaseIx  = model.disc.getDofIx(state, Inf, twoPhase);
            
            ncomp  = model.EOSModel.fluid.getNumberOfComponents();
            cnames = model.EOSModel.fluid.names;
            
            [dx, dy] = deal(zeros(sum(state.nDof), ncomp));
            found = false(ncomp, 1);

            for i = 1:ncomp
                namev = ['v_', cnames{i}];
                namew = ['w_', cnames{i}];
                vix = strcmpi(vars, namev);
                wix = strcmpi(vars, namew);
                if any(wix)
                    found(i) = true;
                    dx(:, i) = rldecode(~pureVapor, state.nDof, 1).*increments{vix};
                    dy(pureVaporIx, i) = increments{vix}(pureVaporIx);
                    dy(twoPhaseIx, i) = increments{wix};
                    [vars, removed] = model.stripVars(vars, {namev, namew});
                    increments = increments(~removed);
                end
            end
            
            n_found = nnz(found);
            if n_found == ncomp-1
                dx(:, ~found) = -sum(dx(:, found), 2);
                dy(:, ~found) = -sum(dy(:, found), 2);
            else
                assert(n_found == 0 || n_found == ncomp);
            end
        end
        
        
        
    end
    
    methods(Access=protected)
        function [ds, dx, vars] = getSatUpdateInternal(model, name, dx, vars, twoPhase, state)
            twoPhaseIx = model.disc.getDofIx(state, Inf, twoPhase);
            dsgix = strcmpi(vars, name);
            if any(dsgix)
                ds_tmp = dx{dsgix};
                if numel(ds_tmp) == numel(twoPhaseIx)
                    ds = zeros(sum(state.nDof), 1);
                    ix = model.disc.getDofIx(state, Inf, twoPhase);
                    ds(ix) = ds_tmp;
                else
                    assert(numel(ds_tmp) == sum(state.nDof));
                    ds = ds_tmp;
                end
                [vars, removed] = model.stripVars(vars, name);
                dx = dx(~removed);
            else
                ds = zeros(sum(state.nDof), 1);
            end
        end
    end
end

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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
