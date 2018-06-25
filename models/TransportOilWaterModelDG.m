classdef TransportOilWaterModelDG < TransportOilWaterModel
    % Two phase oil/water system without dissolution
    properties
        disc
        faceFlux2cellVelocity
    end

    methods
        function model = TransportOilWaterModelDG(G, rock, fluid, varargin)
            
            model = model@TransportOilWaterModel(G, rock, fluid);
            model.disc = [];
            
            model = merge_options(model, varargin{:});
            
            if isempty(model.disc)
                model.disc = DGDiscretization(model, G.griddim);
            end

            
            cellNo    = rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
            faceNo = G.cells.faces(:,1);
            X = G.faces.centroids(faceNo, :) - ...
                G.cells.centroids(cellNo   , :);
            sgn = 1 - 2*(cellNo ~= G.faces.neighbors(faceNo,1));
            
            D1 = sparse(cellNo, faceNo, X(:,1).*sgn, G.cells.num, G.faces.num)./G.cells.volumes;
            D2 = sparse(cellNo, faceNo, X(:,2).*sgn, G.cells.num, G.faces.num)./G.cells.volumes;

            
            model.operators.D = {D1, D2};
            model.operators.faceFlux2cellVelocity = @(v) [D1*v, D2*v];

%             model.degree = 1;
%             model.basis = 'legendre';
% %             model.limiter = 'tvb';
%             
%             model.basis = dgBasis(model.degree, G.griddim, model.basis);
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
            case {'oil', 'sodof'}
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

        function vars = getSaturationVarNames(model)
%             vars = {'sW', 'sO'};
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
            solvedFor = ~strcmpi(saturations, fillsat);
%             ds = zeros(sum(model.disc.nDof), numel(saturations));
            ds = zeros(sum(state.nDof), numel(saturations));
            
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
%             state = model.updateStateFromIncrement(state, ds, problem, 'sdof', inf, inf);
            state = model.updateStateFromIncrement(state, ds, problem, 'sdof', inf, model.dsMaxAbs);
            
            
            
%             state = model.updateStateFromIncrement(state, ds, problem, 'sdof', inf, inf);
            state = model.disc.getCellSaturation(state);
            
            
            if 0
            % Ensure that values are within zero->one interval, and
            % re-normalize if any values were capped
            bad = any((state.s > 1) | (state.s < 0), 2);
            if any(bad)
                cells = find(bad);

%                 ix = model.disc.getDofIx(1, cells);
%                 over  = max(state.s(bad,:)-1,0);
%                 under = min(state.s(bad,:),0);
%                 state.sdof(ix,:) = state.sdof(ix,:) - over - under;
%                 state.sdof(ix,:) = state.sdof(ix,:)./sum(state.sdof(ix,:),2);
                
                s0 = state.s;

                ix = model.disc.getDofIx(1, cells);
%                 state.sdof(ix,:) = min(max(state.s(bad,:),0),1);
                state.sdof(ix,:) = min(max(state.sdof(ix,:),0),1);
                state.sdof(ix,:) = state.sdof(ix,:)./sum(state.sdof(ix,:),2);
                
                ix = model.disc.getDofIx(2:model.disc.basis.nDof, cells);
                state.sdof(ix,:) = 0;
                
                state = model.disc.getCellSaturation(state);
            end
                
            end
            
            if max(model.disc.degree) > 0 && 1
                
                state = model.disc.limiter(state);
                state = model.disc.updateDisc(state);
                
%             sWdof = model.disc.limiter(state.sdof(:,1));
%             sOdof = -sWdof;
%             ix = 1:model.disc.basis.nDof:model.G.cells.num*model.disc.basis.nDof;
%             sOdof(ix) = 1 - sWdof(ix);
% 
%             state.sdof = [sWdof, sOdof];

%                 state = model.disc.getCellSaturation(state);

            end
            
            
            
        end
        
%         function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
%             % Generic update function for reservoir models containing wells.
%             %
%             % SEE ALSO:
%             %   :meth:`ad_core.models.PhysicalModel.updateAfterConvergence`
% 
%             [state, report] = updateAfterConvergence@TransportOilWaterModel(model, state0, state, dt, drivingForces);
%             
%             if model.disc.degree > 0 & 1
% 
%             sWdof = model.disc.limiter(state.sdof(:,1));
%             sOdof = -sWdof;
%             ix = 1:model.disc.basis.nDof:model.G.cells.num*model.disc.basis.nDof;
%             sOdof(ix) = 1 - sWdof(ix);
% 
%             state.sdof = [sWdof, sOdof];
% 
%             state = model.disc.getCellSaturation(state);
% 
%             end
%             
%         end
        
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