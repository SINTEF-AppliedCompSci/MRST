classdef TransportOilWaterModelDG < TransportOilWaterModel
    % Two phase oil/water system without dissolution
    properties
        degree
        basis
        limiter
    end

    methods
        function model = TransportOilWaterModelDG(G, rock, fluid, varargin)
            model = model@TransportOilWaterModel(G, rock, fluid);
            model.degree = 1;
            model.basis = 'legendre';
%             model.limiter = 'tvb';
            
            model = merge_options(model, varargin{:});
            
            model.basis = dgBasis(model.degree, G.griddim, model.basis);
%             model.limiter = dgLimiter(model, model.limiter);

        end

        % --------------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = transportEquationOilWaterDG(state0, state, model,...
                               dt, ...
                               drivingForces,...
                               'solveForOil',   model.conserveOil, ...
                               'solveForWater', model.conserveWater, ...
                               varargin{:});
            
        end
        
        function [fn, index] = getVariableField(model, name)
        % Map variables to state field.
        %
        % SEE ALSO:
        %   :meth:`ad_core.models.PhysicalModel.getVariableField`
        switch(lower(name))
            case {'water', 'swdof'}
                index = 1;
                fn = 'sdof';
            case {'oil', 'swodof'}
                index = 2;
                fn = 'sdof';
            case{'saturation', 'sdof'}
                index = ':';
                fn = 'sdof';
            otherwise
                % This will throw an error for us
                [fn, index] = getVariableField@TransportOilWaterModel(model, name);
        end
        end
        
%         % --------------------------------------------------------------------%
%         function [state, report] = updateState(model, state, problem, dx, drivingForces)
%             % Generic update function for reservoir models containing wells.
%             %
%             % SEE ALSO:
%             %   :meth:`ad_core.models.PhysicalModel.updateState`
% 
%             % Split variables into three categories: Regular/rest variables, saturation
%             % variables (which sum to 1 after updates) and well variables (which live
%             % in wellSol and are in general more messy to work with).
%             [restVars, satVars, wellVars] = model.splitPrimaryVariables(problem.primaryVariables);
% 
%             % Update the wells
%             if isfield(state, 'wellSol')
%                 state.wellSol = model.FacilityModel.updateWellSol(state.wellSol, problem, dx, drivingForces, wellVars);
%             end
% 
%             % Update saturations in one go
%             state  = model.updateSaturations(state, dx, problem, satVars);
% 
%             if ~isempty(restVars)
%                 % Handle pressure seperately
%                 state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMaxRel, model.dpMaxAbs);
%                 state = model.capProperty(state, 'pressure', model.minimumPressure, model.maximumPressure);
%                 restVars = model.stripVars(restVars, 'pressure');
% 
%                 % Update remaining variables (tracers, temperature etc)
%                 for i = 1:numel(restVars)
%                      state = model.updateStateFromIncrement(state, dx, problem, restVars{i});
%                 end
%             end
% 
%             report = [];
%         end

        function vars = getSaturationVarNames(model)
            vars = {'sWdof', 'sOdof'};
            ph = model.getActivePhases();
            vars = vars(ph);
        end
        
        % --------------------------------------------------------------------%
        function state = updateSaturations(model, state, dx, problem, satVars)
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
                [~, satVars] = ...
                    splitPrimaryVariables(model, problem.primaryVariables);
            end
            if isempty(satVars)
                % No saturations passed, nothing to do here.
                return
            end
            % Solution variables should be saturations directly, find the missing
            % link
            saturations = lower(model.getSaturationVarNames);
            fillsat = setdiff(saturations, lower(satVars));
            assert(numel(fillsat) == 1)
            fillsat = fillsat{1};

            % Fill component is whichever saturation is assumed to fill up the rest of
            % the pores. This is done by setting that increment equal to the
            % negation of all others so that sum(s) == 0 at end of update
            nDof = model.basis.nDof;
            solvedFor = ~strcmpi(saturations, fillsat);
            ds = zeros(model.G.cells.num*nDof, numel(saturations));
            
            tmp = 0;
            for i = 1:numel(saturations)
                if solvedFor(i)
                    v = model.getIncrement(dx, problem, saturations{i});
                    ds(:, i) = v;
                    % Saturations added for active variables must be subtracted
                    % from the last phase
                    tmp = tmp - v;
                end
            end
            ds(:, ~solvedFor) = tmp;
            % We update all saturations simultanously, since this does not bias the
            % increment towards one phase in particular.
            state = model.updateStateFromIncrement(state, ds, problem, 'sdof', inf, model.dsMaxAbs);
            
            
            
%             state = model.updateStateFromIncrement(state, ds, problem, 'sdof', inf, inf);
            state = getCellSaturation(model, state);

            if 1
            % Ensure that values are within zero->one interval, and
            % re-normalize if any values were capped
            bad = any((state.s > 1) | (state.s < 0), 2);
            if any(bad)
                ix = (find(bad)-1)*model.basis.nDof + 1;
                over  = max(state.s(bad,:)-1,0);
                under = min(state.s(bad,:),0);
                state.sdof(ix,:) = state.sdof(ix,:) - over - under;
                state.sdof(ix,:) = state.sdof(ix,:)./sum(state.sdof(ix,:),2);
            end
                state = getCellSaturation(model, state);
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