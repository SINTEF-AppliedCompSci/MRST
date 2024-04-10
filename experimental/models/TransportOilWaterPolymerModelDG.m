classdef TransportOilWaterPolymerModelDG < TransportOilWaterModelDG
   
    properties
        polymer
    end
    
    methods
      function model = TransportOilWaterPolymerModelDG(G, rock, fluid, varargin)
            model = model@TransportOilWaterModelDG(G, rock, fluid, varargin{:});
            model.oil     = true;
            model.water   = true;
            model.polymer = true;
            model.conserveWater = false;
            model.conserveOil   = true;
        end

        % ----------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] ...
                = transportEquationOilWaterPolymerDG(state0, state, model, dt, drivingForces, ...
                                  'solveForOil'  , model.conserveOil  , ...
                                  'solveForWater', model.conserveWater, ...
                                  varargin{:}                         );
        end
        % ----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name)
            % Map variables to state field.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.getVariableField`
            switch lower(name)
                case {'cdof'}
                    index = 1;
                    fn = 'cdof';
                case {'c', 'polymer', 'polymermax'}
                    c = model.getComponentNames();
                    index = find(strcmpi(c, 'polymer'));
                    if any(strcmpi(name, {'polymer', 'c'}))
                        fn = 'c';
                    else
                        fn = 'cmax';
                    end
                otherwise
                    [fn, index] = getVariableField@TransportOilWaterModelDG(model, name);
            end
        end
        % ----------------------------------------------------------------%
        function names = getComponentNames(model)
            names = getComponentNames@TransportOilWaterModelDG(model);
            if model.polymer
                names{end+1} = 'polymer';
            end
        end
        % ----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            
            vars = problem.primaryVariables;
            removed = false(size(vars));
            if model.polymer
                % Store the polymer from previous iteration temporarily to
                % use in convergence criteria
                cdof_prev = model.getProp(state, 'polymer');
            end
            
            state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMaxRel, model.dpMaxAbs);
            state = model.capProperty(state, 'pressure', model.minimumPressure, model.maximumPressure);
            [vars, ix] = model.stripVars(vars, 'pressure');
            removed(~removed) = removed(~removed) | ix;
            sWdof = model.getProp(state, 'sWdof');
            dsW   = model.getIncrement(dx, problem, 'sWdof');
            dsO   = -dsW;
                        
            dc = model.getIncrement(dx, problem, 'cDof');
            state = model.updateStateFromIncrement(state, dc, problem, 'cDof');

            nPh = nnz(model.getActivePhases());
            ds = zeros(numel(sWdof), nPh);
            phIndices = model.getPhaseIndices();
            if model.water
                ds(:, phIndices(1)) = dsW;
            end
            if model.oil
                ds(:, phIndices(2)) = dsO;
            end
            state = model.updateStateFromIncrement(state, ds, problem, 'sDof', inf, model.dsMaxAbs);
            if 1
                s = model.disc.getCellMean(state, state.sdof);
                assert(norm(state.s - s) == 0);
                c = model.disc.getCellMean(state, state.cdof);
                assert(norm(state.c - c) == 0);
            end
            
            bad = any((state.s > 1 + model.disc.meanTolerance) ...
                    | (state.s < 0 - model.disc.meanTolerance), 2);
            
            cMax = model.fluid.cmax;
            cBad = any((state.c > cMax*(1 + model.disc.meanTolerance)) ...
                     | (state.c < cMax*(0 - model.disc.meanTolerance)));
            bad = bad | cBad;
                
            if any(bad)
                state.s(bad, :) = min(state.s(bad, :), 1);
                state.s(bad, :) = max(state.s(bad, :), 0);
                state.s(bad, :) = bsxfun(@rdivide, state.s(bad, :), ...
                                              sum(state.s(bad, :), 2));
                state.c(bad)    = min(state.c(bad), cMax);
                state.c(bad)    = max(state.c(bad), 0);
                state = dgLimiter2(model.disc, state, bad);
            end
                
            if model.disc.limitAfterNewtonStep
                [sMin, sMax] = model.disc.getMinMaxSaturation(state);
                bad = sMin < 0 - model.disc.outTolerance | ...
                      sMax > 1 + model.disc.outTolerance;
                if any(bad)
                   state = dgLimiter2(model.disc, state, bad);
                end
                
                if model.disc.jumpTolerance < Inf
                    [cJump, ~, cells] = model.disc.getInterfaceJumps(state.cdof, state);
                    sJump = model.disc.getInterfaceJumps(state.sdof(:,1), state);
                    sJumpTol = model.disc.jumpTolerance;
                    cJumpTol = model.disc.jumpTolerance*model.fluid.cmax;
                    sj = accumarray(cells(:), repmat(sJump,2,1) > sJumpTol) > 0;
                    cj = accumarray(cells(:), repmat(cJump,2,1) > cJumpTol) > 0;
                    j  = sj | cj;
                    jump = false(model.G.cells.num,1);
                    jump(cells(:)) = j(cells(:));
                    bad = jump & state.degree > 0;
                    if any(bad)
                        state = dgLimiter2(model.disc, state, bad);
                    end
                end
            end
            
            %  We have explicitly dealt with rs/rv properties, remove from list
            %  meant for autoupdate.
            [vars, ix] = model.stripVars(vars, {'sWdof', 'sOdof', 'cDof'});
            removed(~removed) = removed(~removed) | ix;

            % We may have solved for a bunch of variables already if we had
            % disgas / vapoil enabled, so we remove these from the
            % increment and the linearized problem before passing them onto
            % the generic reservoir update function.
            problem.primaryVariables = vars;
            dx(removed) = [];

            % Parent class handles almost everything for us
            [state, report] = updateState@TransportOilWaterModelDG(model, state, problem, dx, drivingForces);
            
            if model.polymer && 0
                
                state = model.capProperty(state, 'c', 0, model.fluid.cmax);
                state.c_prev = cdof_prev;

                % Shear Thinning Report
                % We (may) have stored the shear thinning report
                % temporarily in the state structure. We move this over to
                % the report structure instead. The reason for this is that
                % there is no report returned from the equations.
                if isfield(state, 'ShearThinningReport')
                    report.ShearThinning = state.ShearThinningReport;
                    state = rmfield(state, 'ShearThinningReport');
                end
            end
            
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@TransportOilWaterModelDG(model, state0, state, dt, drivingForces);
            if model.polymer
%                 state = model.capProperty(state, 'c', 0, model.fluid.cmax);
                c     = model.getProp(state, 'polymer');
                cmax  = model.getProp(state, 'polymermax');
                state = model.setProp(state, 'polymermax', max(cmax, c));

                if isfield(state, 'c_prev')
                    % Remove the temporary field used for convergence
                    state = rmfield(state, 'c_prev');
                end
            end
        end

    end
    
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
