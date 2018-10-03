classdef TransportBlackOilModelDG < TransportBlackOilModel
    % Two phase oil/water system without dissolution with discontinuous
    % Galerking discretization
    
    properties
        disc % DG discretization
    end

    methods
        function model = TransportBlackOilModelDG(G, rock, fluid, varargin)
            
            model = model@TransportBlackOilModel(G, rock, fluid);
            model.disc = [];
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
                case {'water', 'swdof'}
                    index = 1;
                    fn = 'sdof';
                case {'oil', 'sodof'}
                    index = 2;
                    fn = 'sdof';
                case {'gas', 'sgdof'}
                    index = 3;
                    fn = 'sdof';
                case{'saturation', 'sdof'}
                    index = ':';
                    fn = 'sdof';
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@TransportBlackOilModel(model, name);
            end
        end

        % ----------------------------------------------------------------%
        function vars = getSaturationVarNames(model)
            vars = {'sWdof', 'sOdof', 'sGdof'};
            ph = model.getActivePhases();
            vars = vars(ph);
        end
        
        %-----------------------------------------------------------------%
        function integrand = cellIntegrand(model, fun, x, cellNo, state, state0, sdof, sdof0, f)            
            
            % Evaluate saturations and fractional flow at cubature points
            s  = model.disc.evaluateSaturation(x, cellNo, sdof , state );
            s0 = model.disc.evaluateSaturation(x, cellNo, sdof0, state0);
            f = f(s, 1-s, cellNo, cellNo);
            integrand = @(psi, grad_psi) fun(s, s0, f, cellNo, psi, grad_psi);
            
        end
        
        %-----------------------------------------------------------------%
        function integrand = faceIntegrand(model, fun, x, faceNo, cellNo, T, vT, g, mob, state, sdof, f)
            % TODO: See if upwind direction changes a lot during nonlinear
            % solution. If so, maybe fix upstream after iteratio n as in "A
            % fully-coupled upwind discontinuous ..."
            
            % Get upstream cells and upstram saturations at given
            % quadrature points
            [flag_v, flag_G, upCells_v, upCells_G, s_v, s_G] ...
                = model.disc.getSaturationUpwind(faceNo, x, T, vT, g, mob, sdof, state);
            
            f_v = f(s_v{1}, 1-s_v{2}, upCells_v{1}, upCells_v{2});
            f_G = f(s_G{1}, 1-s_G{2}, upCells_G{1}, upCells_G{2});
            
            integrand = @(psi) fun(s_v{1}, s_G{1}, f_v, f_G, cellNo, upCells_v{1}, upCells_G{1}, faceNo, psi);

        end
        
        %-----------------------------------------------------------------%
        function integrand = faceIntegrandBC(model, fun, x, faceNo, cellNo, bc, isInj, globFace2BCface, state, sdof, f)
            
            % Get upstream cells and upstram saturations at given
            % quadrature points
            
            
            locFaceNo = globFace2BCface(faceNo);
            
            sL = bc.sat(:,1);
            sL = sL(locFaceNo);
            sR = model.disc.evaluateSaturation(x, cellNo, sdof, state);
            s = sL.*isInj(locFaceNo) + sR.*(~isInj(locFaceNo));
            f = f(s, 1-s, cellNo, cellNo);
            
            integrand = @(psi) fun(s, f, cellNo, faceNo, psi);
            
            
        end
        
        % ----------------------------------------------------------------%
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
            active = ~model.G.cells.ghost;
            ix = model.disc.getDofIx(state, Inf, active);
            for i = 1:numel(saturations)
                if solvedFor(i)
                    v = model.getIncrement(dx, problem, saturations{i});
                    ds(ix, i) = v;
                    % Saturations added for active variables must be subtracted
                    % from the last phase
                    tmp = tmp - v;
                end
            end
            ds(ix, ~solvedFor) = tmp;
            % We update all saturations simultanously, since this does not bias the
            % increment towards one phase in particular.
%             state = model.updateStateFromIncrement(state, ds, problem, 'sdof', inf, inf);
            state = model.updateStateFromIncrement(state, ds, problem, 'sdof', inf, model.dsMaxAbs);
            
            
            
%             state = model.updateStateFromIncrement(state, ds, problem, 'sdof', inf, inf);
            state.s = model.disc.getCellSaturation(state);
            
            if 1
                state = model.disc.limiter(state);
            end
            
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            % Generic update function for reservoir models containing wells.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.updateAfterConvergence`

            [state, report] = updateAfterConvergence@TransportBlackOilModel(model, state0, state, dt, drivingForces);
            
            if model.disc.degree > 0 && model.disc.jumpTolerance < Inf
                % Cells with interface jumps larger than threshold
                [jumpVal, ~, cells] = model.disc.getInterfaceJumps(state.sdof(:,1), state);
                j = accumarray(cells(:), repmat(jumpVal,2,1) > model.disc.jumpTolerance) > 0;
                jump = false(model.G.cells.num,1);
                jump(cells(:))          = j(cells(:));
                jump(state.degree == 0) = false;
                if any(jump)
                    state = dgLimiter(model.disc, state, jump, 'tvb');
                end
%                 state = model.disc.limiter(state);
            end
            
            if model.disc.degree > 0 && 1
                state = dgLimiter(model.disc, state, true(model.G.cells.num,1), 'scale');
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