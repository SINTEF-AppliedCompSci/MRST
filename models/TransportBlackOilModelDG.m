classdef TransportBlackOilModelDG < TransportBlackOilModel
    % Two phase oil/water system without dissolution with discontinuous
    % Galerking discretization
    
    properties
        disc % DG discretization
        tryMaxDegree
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
            model = merge_options(model, varargin{:});
            
            % Construct discretization
            if isempty(model.disc)
                model.disc = DGDiscretization(model, G.griddim);
            end

        end

        % ----------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] ...
                = transportEquationOilWaterDG(state0, state, model, dt, drivingForces, ...
                                  'solveForOil'  , model.conserveOil  , ...
                                  'solveForWater', model.conserveWater, ...
                                  'solveForGas'  , model.conserveGas  , ...
                                  varargin{:}                         );
            
        end
        
        % ----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name)
            % Map variables to state field.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.getVariableField`
            switch(lower(name))
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
        
        % ----------------------------------------------------------------%
        function state = updateSaturations(model, state, dx, problem, satDofVars)
            % Update of phase-saturations
            %
            % SYNOPSIS:
            %   state = model.updateSaturations(state, dx, problem, satVars)
            %
            % DESCRIPTION:
            %   Update saturations (likely state.s) under the constraint that
            %   the sum of volume fractions is always equal to 1. This
            %   assumes that we have solved for n - 1 phases when n phases
            %   are present.
            %
            % PARAMETERS:
            %   model   - Class instance
            %   state   - State to be updated
            %   dx      - Cell array of increments, some of which correspond 
            %             to saturations
            %   problem - `LinearizedProblemAD` class instance from which `dx`
            %             was obtained.
            %   satVars - Cell array with the names of the saturation
            %             variables.
            %
            % RETURNS:
            %   state - Updated state with saturations within physical
            %           constraints.
            %
            % SEE ALSO:
            %   `splitPrimaryVariables`

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
                % Fill component is whichever saturation is assumed to fill up the rest of
                % the pores. This is done by setting that increment equal to the
                % negation of all others so that sum(s) == 0 at end of update
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
            state   = model.updateStateFromIncrement(state, ds, problem, 'sdof', inf, model.dsMaxAbs);
            state.s = model.disc.getCellSaturation(state);
            
            if nFill == 1
                bad = any((state.s > 1 + model.disc.outTolerance) ...
                        | (state.s < 0 - model.disc.outTolerance), 2);
                if any(bad)
                    state.s(bad, :) = min(state.s(bad, :), 1);
                    state.s(bad, :) = max(state.s(bad, :), 0);
                    state.s(bad, :) = bsxfun(@rdivide, state.s(bad, :), ...
                                                  sum(state.s(bad, :), 2));
                    state = dgLimiter(model.disc, state, bad, 'kill');
                end
            else
                bad = any(state.s < 0 - model.disc.outTolerance, 2);
                 if any(bad)
                    state.s(bad, :) = max(state.s(bad, :), 0);
                    state = dgLimiter(model.disc, state, bad, 'kill');
                 end
            end

            if model.disc.limitAfterNewtonStep
                % Limit solution
                state = model.disc.limiter(model, state, [], true);
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