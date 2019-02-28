classdef TransportBlackOilModelDG < TransportBlackOilModel
    % Two phase oil/water system without dissolution with discontinuous
    % Galerking discretization
    
    properties
        disc         % DG discretization
        tryMaxDegree % 
    end

    methods
        % ----------------------------------------------------------------%
        function model = TransportBlackOilModelDG(G, rock, fluid, varargin)
            
            model = model@TransportBlackOilModel(G, rock, fluid);
            model.disc = [];
            model.tryMaxDegree = true;
            % If we use reordering, this tells us which cells are actually
            % part of the discretization, and which cells that are included
            % to get fluxes correct
            model.G.cells.ghost = false(G.cells.num,1);
            [model, discArgs] = merge_options(model, varargin{:});
            
            % Construct discretization
            if isempty(model.disc)
                model.disc = DGDiscretization(model, discArgs{:});
            end

        end

        % ----------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] ...
                = transportEquationBlackOilDG(state0, state, model, dt, drivingForces, ...
                                  'solveForOil'  , model.conserveOil  , ...
                                  'solveForWater', model.conserveWater, ...
                                  'solveForGas'  , model.conserveGas  , ...
                                  varargin{:}                         );
            
        end
        
        % ----------------------------------------------------------------%
        function kr = evaluateRelPerm(model, phase, sat, varargin)
            
            switch phase
                case 'W'
                    name = 'W';
                case 'O'
                    name = ['O_', model.getPhaseNames()];
                case 'G'
                    name = 'G';
            end
            
            fn = ['relPerm', name];
            kr = model.(fn)(sat, varargin{:});
            
        end
        
        % ----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name)
            % Map variables to state field.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.getVariableField`
            switch lower(name)
                case {'swdof'}
                    index = 1;
                    fn = 'sdof';
                case {'sodof'}
                    index = 2;
                    fn = 'sdof';
                case {'sgdof'}
                    index = 3;
                    fn = 'sdof';
                case{'sdof'}
                    index = ':';
                    fn = 'sdof';
                case {'rsdof', 'rvdof'}
                    fn = lower(name);
                    index = 1;
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@TransportBlackOilModel(model, name);
            end
        end

        % ----------------------------------------------------------------%
        function vars = getDGDofVarNames(model)
            vars = {'swdof', 'sodof', 'sgdof'};
            ph = model.getActivePhases();
            vars = vars(ph);
        end
        
        % ----------------------------------------------------------------%
        function [restVars, satDofVars, wellVars] = splitPrimaryVariables(model, vars)
            % Split cell array of primary variables into grouping
            % SYNOPSIS:
            %   [restVars, satVars, wellVars] = model.splitPrimaryVariables(vars)
            %
            % DESCRIPTION:
            %   Split a set of primary variables into three groups:
            %   Well variables, saturation variables and the rest. This is
            %   useful because the saturation variables usually are updated
            %   together, and the well variables are a special case.
            %
            % PARAMETERS:
            %   model - Class instance.
            %   vars  - Cell array with names of primary variables
            %
            % RETURNS:
            %   restVars - Names of variables that are not saturations or
            %              belong to the wells.
            %   satVars  - Names of the saturation variables present in `vars`.
            %   wellVars - Names of the well variables present in `vars`
            
            wellvars = model.FacilityModel.getPrimaryVariableNames();
            isSatDof = cellfun(@(x) any(strcmpi(model.getDGDofVarNames, x)), vars);
            isWells  = cellfun(@(x) any(strcmpi(wellvars, x)), vars);

            wellVars   = vars(isWells);
            satDofVars = vars(isSatDof);

            restVars = vars(~isSatDof & ~isWells);
        end
        
        % --------------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            vars = problem.primaryVariables;
            removed = false(size(vars));
            if model.disgas || model.vapoil
                % The VO model is a bit complicated, handle this part
                % explicitly.
                state0 = state;

                state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMaxRel, model.dpMaxAbs);
                state = model.capProperty(state, 'pressure', model.minimumPressure, model.maximumPressure);

                [vars, ix] = model.stripVars(vars, 'pressure');
                removed(~removed) = removed(~removed) | ix;

                % Black oil with dissolution
                sOdof = model.getProp(state, 'sOdof');
                sGdof = model.getProp(state, 'sGdof');
                if model.water
                    sWdof = model.getProp(state, 'sWdof');
                    dsW = model.getIncrement(dx, problem, 'sWdof');
                else
                    sWdof = 0;
                    dsW = 0;
                end
                % Magic status flag, see inside for doc
                [sW, sO, sG] = model.disc.getCellMean(state, sWdof, sOdof, sGdof);
                st = model.getCellStatusVO(state0, sO, sW, sG);
                stt = st;
                st  = cellfun(@(st) expandSingleValue(st, model.G), st , 'unif', false);
                st  = cellfun(@(st) rldecode(st, state.nDof, 1), st, 'unif', false);
                
                dr = model.getIncrement(dx, problem, 'xDof');
                % Interpretation of "gas" phase varies from cell to cell, remove
                % everything that isn't sG updates
                dsG = st{3}.*dr - st{2}.*dsW;

                if model.disgas
                    state = model.updateStateFromIncrement(state, st{1}.*dr, problem, ...
                                                           'rSdof', model.drsMaxRel, model.drsMaxAbs);
                    state.rs = model.disc.getCellMean(state, state.rsdof);
                end

                if model.vapoil
                    state = model.updateStateFromIncrement(state, st{2}.*dr, problem, ...
                                                           'rVdof', model.drsMaxRel, model.drsMaxAbs);
                    state.rv = model.disc.getCellMean(state, state.rvdof);
                end

                dsO = -(dsG + dsW);
                nPh = nnz(model.getActivePhases());

                ds = zeros(numel(sOdof), nPh);
                phIndices = model.getPhaseIndices();
                if model.water
                    ds(:, phIndices(1)) = dsW;
                end
                if model.oil
                    ds(:, phIndices(2)) = dsO;
                end
                if model.gas
                    ds(:, phIndices(3)) = dsG;
                end

                state = model.updateStateFromIncrement(state, ds, problem, 'sDof', inf, model.dsMaxAbs);
                state.s = model.disc.getCellSaturation(state);
                % We should *NOT* be solving for oil saturation for this to make sense
                assert(~any(strcmpi(vars, 'sOdof')));
                state = computeFlashBlackOilDG(state, state0, model, stt);
                
                bad = any((state.s > 1 + model.disc.meanTolerance) ...
                        | (state.s < 0 - model.disc.meanTolerance), 2);    
                    
                if any(bad)
                    state.s(bad, :) = min(state.s(bad, :), 1);
                    state.s(bad, :) = max(state.s(bad, :), 0);
                    state.s(bad, :) = bsxfun(@rdivide, state.s(bad, :), ...
                                                  sum(state.s(bad, :), 2));
                    state = dgLimiter2(model.disc, state, bad);
                end
                sT = sum(state.s, 2);
                state.sdof = bsxfun(@rdivide, state.sdof, rldecode(sT, state.nDof,1));
                state.s    = bsxfun(@rdivide, state.s, sT);
                
                if model.disc.limitAfterNewtonStep
                    [sMin, sMax] = model.disc.getMinMaxSaturation(state);
                    bad = sMin < 0 - model.disc.outTolerance | ...
                          sMax > 1 + model.disc.outTolerance;
                    if any(bad)
                       state = dgLimiter2(model.disc, state, bad);
                    end
                end
                %  We have explicitly dealt with rs/rv properties, remove from list
                %  meant for autoupdate.
                [vars, ix] = model.stripVars(vars, {'sWdof', 'sOdof', 'sGdof', 'rSdof', 'rVdof', 'xDof'});
                removed(~removed) = removed(~removed) | ix;

            end

            % We may have solved for a bunch of variables already if we had
            % disgas / vapoil enabled, so we remove these from the
            % increment and the linearized problem before passing them onto
            % the generic reservoir update function.
            problem.primaryVariables = vars;
            dx(removed) = [];

            % Parent class handles almost everything for us
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
        end
        
        
        % ----------------------------------------------------------------%
        function state = updateSaturations(model, state, dx, problem, satDofVars)

            if nargin < 5
                % Get the saturation names directly from the problem
                [~, satDofVars] = ...
                    splitPrimaryVariables(model, problem.primaryVariables);
            end
            if isempty(satDofVars)
                % No saturations passed, nothing to do here.
                return
            end
            % Solution variables should be saturations directly, find the missing
            % link
            saturations = lower(model.getDGDofVarNames);
            
            fillsat = setdiff(saturations, lower(satDofVars));
            nFill = numel(fillsat);
            assert(nFill == 0 || nFill == 1)
            if nFill == 1
                % Fill component is whichever saturation is assumed to fill
                % up the rest of the pores. This is done by setting that
                % increment equal to the negation of all others so that
                % sum(s) == 0 at end of update
                fillsat = fillsat{1};
                solvedFor = ~strcmpi(saturations, fillsat);
            else
                % All saturations are primary variables. Sum of saturations is
                % assumed to be enforced from the equation setup
                solvedFor = true(numel(saturations), 1);
            end
            ds = zeros(sum(state.nDof), numel(saturations));
            
            tmp = 0;
            active = ~model.G.cells.ghost;
            ix = model.disc.getDofIx(state, Inf, active);
            for phNo = 1:numel(saturations)
                if solvedFor(phNo)
                    v = model.getIncrement(dx, problem, saturations{phNo});
                    ds(ix, phNo) = v;
                    if nFill > 0
                        % Saturations added for active variables must be subtracted
                        % from the last phase
                        tmp = tmp - v;
                    end
                end
            end
            ds(ix, ~solvedFor) = tmp;
            % We update all saturations simultanously, since this does not bias the
            % increment towards one phase in particular.
            state   = model.updateStateFromIncrement(state, ds, problem, 'sdof', Inf, model.dsMaxAbs);
            state.s = model.disc.getCellSaturation(state);

%             ix = any(abs(ds)>model.dsMaxAbs,2);
%             alph = model.dsMaxAbs./max(abs(ds(ix,:)), [], 2);
%             
%             ds(ix,:) = alph.*ds(ix,:);
%             state   = model.updateStateFromIncrement(state, ds, problem, 'sdof', Inf, Inf);
%             state.s = model.disc.getCellSaturation(state);            
            
            if nFill == 1
                
                if 1
                bad = any((state.s > 1 + model.disc.meanTolerance) ...
                        | (state.s < 0 - model.disc.meanTolerance), 2);
                
                else
                [smin, smax] = model.disc.getMinMaxSaturation(state);
                over  = smax > 1 + model.disc.meanTolerance;
                under = smin < 0 - model.disc.meanTolerance;
                bad = over | under;
                end
                    
                    
                if any(bad)
                    state.s(bad, :) = min(state.s(bad, :), 1);
                    state.s(bad, :) = max(state.s(bad, :), 0);
                    state.s(bad, :) = bsxfun(@rdivide, state.s(bad, :), ...
                                                  sum(state.s(bad, :), 2));
                    state = dgLimiter(model.disc, state, bad, 's', 'kill');
                end
            else
                bad = any(state.s < 0 - model.disc.meanTolerance, 2);
                 if any(bad)
                    state.s(bad, :) = max(state.s(bad, :), 0);
                    state = dgLimiter(model.disc, state, bad, 's', 'kill');
                 end
            end

            if model.disc.limitAfterNewtonStep
                % Limit solution
                state = model.disc.limiter(model, state, [], true);
            end
            
        end
        
        function [state, val, val0] = updateStateFromIncrement(model, state, dx, problem, name, relchangemax, abschangemax)
            
            if strcmpi(name, 'sdof') && abschangemax < Inf && 1
                
                if 1
                    s0 = state.s;

                    st = state;
                    st.sdof = st.sdof + dx;

                    s  = model.disc.getCellSaturation(st);
                    ds = abs(s - s0);
                    outside = any(ds > abschangemax,2);
                    alpha = abschangemax./max(ds(outside), [], 2);
    %                 ix = rldecode(find(outside), state.nDof(outside), 1);
                    dofIx = model.disc.getDofIx(state, Inf, outside);
                    dx(dofIx,:) = rldecode(alpha, state.nDof(outside), 1).*dx(dofIx,:);

                    sdof = state.sdof + dx;
                    state = model.setProp(state, name, sdof);
                    
                elseif 0
                    
                    
                
                else
                    
                    dxOutside = zeros(model.G.cells.num,1);
                    for dofNo = 1:model.disc.basis.nDof
                        ix = model.disc.getDofIx(state, dofNo, Inf, true);
                        
                        outside = false(model.G.cells.num,1);
                        outside(ix>0) = any(abs(dx(ix(ix>0),:)) > abschangemax,2);

                        dxOutside(outside) = max(abs(dx(ix(outside))), [], 2);
                    end
                    
                    outside = dxOutside > 0;
                    alpha = ones(model.G.cells.num,1);
                    alpha(outside) = abschangemax./dxOutside(outside);
                    dx = dx.*rldecode(alpha, state.nDof, 1);
                    sdof = state.sdof + dx;
                    state = model.setProp(state, name, sdof);
                
                end

            else
                [state, val, val0] ...
                    = updateStateFromIncrement@TransportBlackOilModel(model, state, dx, problem, name, relchangemax, abschangemax);
            end
            
        end
        

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            % Generic update function for reservoir models containing wells.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.updateAfterConvergence`

            % Let base class do what it want
            [state, report] = updateAfterConvergence@TransportBlackOilModel(model, state0, state, dt, drivingForces);
            if model.disc.limitAfterConvergence
                % Postprocess using limiter(s)
                state = model.disc.limiter(model, state, state0, false);    
            end
        end
        
        function krW = relPermW(model, sat, varargin)
            sW  = sat{model.getPhaseIndex('W')};
            krW = model.fluid.krW(sW, varargin{:});
        end

        function krO = relPermO_O(model, sat, varargin)
            sO  = sat{model.getPhaseIndex('O')};
            krO = model.fluid.krO(sO, varargin{:});
        end

        function krO = relPermO_WO(model, sat, varargin)
            sO = sat{model.getPhaseIndex('O')};
            if isfield(model.fluid, 'krO')
                krO = model.fluid.krO(sO, varargin{:});
            else
                krO = model.fluid.krOW(sO, varargin{:});
            end
        end

        function krO = relPermO_OG(model, sat, varargin)
            sO = sat{model.getPhaseIndex('O')};
            if isfield(model.fluid, 'krO')
                krO = model.fluid.krO(sO, varargin{:});
            else
                krO = model.fluid.krOG(sO, varargin{:});
            end
        end

        function krO = relPermO_WOG(model, sat, varargin)
            sWcon = 0;
            fluid = model.fluid;
            if isfield(fluid, 'sWcon')
                if isempty(varargin) || numel(fluid.sWcon) == 1
                    sWcon = fluid.sWcon;
                else
                    assert(strcmp(varargin{1}, 'cellInx'))
                    sWcon = fluid.sWcon(varargin{2});
                end
            end
            sW = sat{model.getPhaseIndex('W')};
            sO = sat{model.getPhaseIndex('O')};
            sG = sat{model.getPhaseIndex('G')};
            sWcon = min(sWcon, double(sW)-1e-5);

            wW = (sW-sWcon)./(sG+sW-sWcon);
            wG = 1-wW;

            krOW = fluid.krOW(sO, varargin{:});
            krOG = fluid.krOG(sO,  varargin{:});
            krO  = wG.*krOG + wW.*krOW;

        end      

        function krG = relPermG(model, sat, varargin)
            sG  = sat{model.getPhaseIndex('G')};
            krG = model.fluid.krG(sG, varargin{:});
        end
        
    end
end

function v = expandSingleValue(v,G)
    if numel(double(v)) == 1
        v = v*ones(G.cells.num,1);
    end
end

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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