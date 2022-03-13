classdef CompositionalMSPressureModel < MultiscalePressureModel
    % Compositional multiscale model (which also works for regular
    % sequential without compositional effects)
    properties
        resetBasisAfterConvergence
        regularizeReconstructionUpdate
        reconstructWellFlux
        reduceWellSystem
        useGetEquationsForFlux = true;
    end
    
    methods
        
        function model = CompositionalMSPressureModel(G, rock, fluid, pmodel, mssolver, varargin)
            model = model@MultiscalePressureModel(G, rock, fluid, pmodel, mssolver);
            model.resetBasisAfterConvergence = false;
            model.regularizeReconstructionUpdate = false;
            assert(~isa(pmodel, 'CompositionalMSPressureModel'));
            model.useGetEquationsForFlux = ~isa(pmodel, 'PressureModel');
            model.dpMaxRel = inf;
            model.reconstructWellFlux = false;
            model.reduceWellSystem = false;
            model = merge_options(model, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = model.pressureModel.getEquations(state0, state, dt, drivingForces, varargin{:}, 'resOnly', false);
            problem = problem.assembleSystem();
            conv = model.checkConvergence(problem);
            if ~all(conv)
                state.tmp.sys = problem;
            end
        end

        function [convergence, values, names] = checkConvergence(model, problem, varargin)
            [convergence, values, names] = checkConvergence(model.pressureModel, problem, varargin{:});
            if model.coarseScaleConvergence
                sub = strcmpi(names, 'pressure');
                if any(sub)
                    if isempty(model.multiscaleSolver.basis)
                        R = double(controlVolumeRestriction(model.multiscaleSolver.coarsegrid.partition));
                    else
                        R = model.multiscaleSolver.basis.R;
                    end
                    values(sub) = norm(R*problem.b(1:model.G.cells.num), inf);
                    convergence(sub) = values(sub) < model.pressureModel.nonlinearTolerance;
                    names{1} = 'Pressure';
                end
                sub = strcmpi(problem.equationNames, 'volclosure');
                twoPhase = problem.state.flag == 0;
                if any(sub) && any(twoPhase)
                    part = model.multiscaleSolver.coarsegrid.partition;
                    
                    f_eq = double(problem.equations{sub});
                    f_eq = f_eq.*model.operators.pv(twoPhase);
                    
                    f_eq_a = zeros(size(model.operators.pv));
                    f_eq_a(twoPhase) = f_eq;
                    
                    pv = model.operators.pv.*twoPhase;
                    values(sub) = norm(accumarray(part, f_eq_a)./accumarray(part, pv), inf);
                    convergence(sub) = values(sub) < model.pressureModel.nonlinearTolerance;
                end
            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, forces)
            [state, report] = updateAfterConvergence@MultiscalePressureModel(model, state0, state, dt, forces);
            if isfield(state, 'tmp')
                state = rmfield(state, 'tmp');
            end
            if model.resetBasisAfterConvergence
                model.multiscaleSolver.basis = [];
            end
        end
    
        function [state, varargout] = updateState(model, state, problem, dx, varargin)
            varargout = cell(1, nargout-1);
            p0 = state.pressure;
            isPressure = strcmpi(problem.primaryVariables, 'pressure');
            cg = model.multiscaleSolver.coarsegrid;
            
            dp = dx{isPressure};
            if isfinite(model.dpMaxRel)
                [~, w] = model.limitUpdateRelative(dp, state.pressure, model.dpMaxRel);
                w = min(w);
                for i = 1:numel(dx)
                    dx{i} = dx{i}*w;
                end
            end
            if model.multiscaleSolver.maxIterations == 0
                dx0 = dx;
                B = model.multiscaleSolver.basis.B;
%                 centers = cg.cells.centers;
                [~, centers] = max(B);
                p = state.pressure + dp;
                p_c = p(centers);
                p_ms = B*p_c;
                dp_ms = p_ms - p0;
                if 1
                    dx{isPressure} = dp_ms;
                else
                    [A, b, B, C, D, E, f, h] = reduceSystem(problem.A, problem.b, numel(p));
                    s = E\(h - D*dp_ms);
                    dp_ms = [dp_ms; s];
                    increments = problem.processResultAfterSolve(dp_ms, []);
                    sol = LinearSolverAD();
                    dx = sol.storeIncrements(problem, increments);
                end
                model.pressureModel.dpMaxRel = inf;
                model.pressureModel.dpMaxAbs = inf;
            end
            [state, varargout{:}] = model.pressureModel.updateState(state, problem, dx, varargin{:});            
            state.tmp.dpIncrement = state.pressure - p0;
            state.tmp.dpRaw = dp;
            if model.coarseScaleConvergence
                if isfield(state, 'dpRel')
                    state.dpRel = state.dpRel(cg.cells.centers);
                end
            end
        end

        % Utility functions 
        function [state, report] = reconstructFluxes(model, state0, stateConverged, dt, forces)
            CG = model.multiscaleSolver.coarsegrid;
            nc = CG.parent.cells.num;
            
            timer = tic();
            % Property pressure is the converged pressure carried forward
            % to transport
            propPressure = stateConverged.pressure;
            % We are a little bit strict on defining the last pressure
            % increment since the nonlinear solver and update logic may
            % modify the increment in ways that are incompatible with the
            % linearized coarse scale finite-volume scheme. We re-caculate
            % this value through the stored increments.
            dpMS = stateConverged.tmp.dpRaw;
            p0 = stateConverged.pressure - stateConverged.tmp.dpIncrement;
            if 1
                % This was the last linear system used to solve a step.
                problem = stateConverged.tmp.sys;
                A = problem.A;
                b = problem.b;
            else
                extra = {};
                problem = model.pressureModel.getEquations(state0, stateConverged, dt, forces,...
                    'iteration', inf, 'staticwells', true, extra{:});
                if isa(problem, 'PressureReducedLinearSystem')
                    A = problem.A;
                    b = problem.b;
                else
                    A = problem.equations{1}.jac{1};
                    b = -problem.equations{1}.val;
                end
            end
            if model.reduceWellSystem
                solver = model.multiscaleSolver;
                [A, b, lsys] = solver.reduceLinearSystem(A, b);
            else
                A = A(1:nc, 1:nc);
                b = b(1:nc);
            end

            [stateCoarse, stateFlux] = deal(stateConverged);
            % Pressure at second to last multiscale solution, plus the raw
            % increment produced by the multiscale solver.
            stateCoarse.pressure = p0 + dpMS;
            if model.regularizeReconstructionUpdate
                [dpReconstruct, t_solve] = reconstructPressureNormalized(CG, dpMS, A, b);
            else
                [dpReconstruct, t_solve] = reconstructPressure(CG.partition, dpMS, A, b, true);
            end
            % Update the internal domain pressure with the reconstructed
            % increment.
            if model.reduceWellSystem
                pmodel = model.pressureModel;
                pmodel.dpMaxRel = inf;
                pmodel.dpMaxAbs = inf;
                dpReconstruct = solver.recoverLinearSystem(dpReconstruct, lsys);
                dx = solver.storeIncrements(problem, dpReconstruct);
                [stateFlux, report] = pmodel.updateState(stateFlux, problem, dx, forces);
            else
                stateFlux.pressure =  p0 + dpReconstruct(1:nc);
            end
            % Store fluxes and use it the corect places
            stateFlux   = setFluxes(model, state0, stateFlux,   dt, forces, propPressure);
            stateCoarse = setFluxes(model, state0, stateCoarse, dt, forces, propPressure);

            flux = stateFlux.flux;
            flux(CG.faces.fconn, :) = stateCoarse.flux(CG.faces.fconn, :);
            state = stateConverged;
            state.flux = flux;
            if model.reconstructWellFlux
                for i = 1:numel(state.wellSol)
                    state.wellSol(i).flux = stateFlux.wellSol(i).flux;
                end
            end            
            % Use property pressure for transport
            state.pressure = propPressure;

            report.reconstructionTime = toc(timer);
            report.reconstructionSolver = t_solve;
        end
        function state = setFluxes(model, state0, state, dt, forces, propsPressure)
            m = model.pressureModel;
            if model.useGetEquationsForFlux
                [~, state] = m.getEquations(state0, state, dt, forces, ...
                    'propsPressure', propsPressure, 'resOnly', true, 'iteration', inf);
            else
                p_flux = state.pressure;
                evalstate = m.initStateFunctionContainers(state);
                m.getProp(evalstate, 'PressureGradient'); % Trigger evaluation
                evalstate.pressure = propsPressure;
                f = model.pressureModel.getProp(evalstate, 'PhaseFlux');
                nph = numel(f);
                state.flux = zeros(model.G.faces.num, nph);
                state.flux(model.operators.internalConn, :) = value(f);

                if ~isempty(forces.W) && model.reconstructWellFlux
                    fm = m.FacilityModel;
                    [mob] = m.getProps(evalstate, 'Mobility');
                    p = p_flux;
                    tmp = fm.FacilityFluxDiscretization.getStateFunction('FacilityWellMapping');
                    map = tmp.evaluateOnDomain(fm, state);
                    cdp = vertcat(state.wellSol.cdp);
                    pw = vertcat(state.wellSol(map.active).bhp);

                    dp = pw(map.perf2well) + cdp - p(map.cells);

                    T = vertcat(map.W.WI);
                    mobt = value(mob);
                    qw = T.*dp.*sum(mobt(map.cells, :), 2);
                    for i = 1:numel(map.active)
                        ix = map.active(i);
                        state.wellSol(ix).flux = qw(map.perf2well == ix);
                    end
                end
            end
        end
    end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
