classdef MultiscalePressureModel < ReservoirModel
    properties
        pressureModel
        multiscaleSolver
        blockIndices
        useMex
        useReconstructedPressure
        coarseScaleConvergence
        reconstruct
    end
    
    methods
        
        function model = MultiscalePressureModel(G, rock, fluid, pmodel, mssolver, varargin)
            model = model@ReservoirModel(G, rock, fluid);
            [model.water, model.oil, model.gas] = deal(pmodel.water, pmodel.oil, pmodel.gas);
            model.useMex = false;
            model.pressureModel = pmodel;
            model.multiscaleSolver = mssolver;
            model.useReconstructedPressure = false;
            model.coarseScaleConvergence = false;
            model.reconstruct = true;
            
            model = merge_options(model, varargin{:});
            
            if model.useMex
                CG = mssolver.coarsegrid;
                mrstModule add mpfa
                tmp = accumarray(CG.partition, ones(G.cells.num, 1));
                [i, j] = blockDiagIndex(tmp);

                model.blockIndices = struct('I', i, 'J', j, 'blocksz', tmp);
            end
        end
        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolve, nonlinsolve, iteration, varargin)
%             [state, report] = model.pressureModel.stepFunction(state, state0, dt, drivingForces, model.multiscaleSolver, nonlinsolve, iteration, varargin{:});
            [state, report] = stepFunction@PhysicalModel(model, state, state0, dt, drivingForces, model.multiscaleSolver, nonlinsolve, iteration, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = model.pressureModel.getEquations(state0, state, dt, drivingForces, varargin{:});
        end
        
        function [vars, names, origin] = getPrimaryVariables(model, state)
            [vars, names, origin] = model.pressureModel.getPrimaryVariables(state);
        end

        function state = initStateAD(model, state, vars, names, origin)
            state = model.pressureModel.initStateAD(state, vars, names, origin);
        end
        
        function [eqs, names, types, state] = getModelEquations(pmodel, state0, state, dt, drivingForces)
            [eqs, names, types, state] = model.pressureModel.getModelEquations(state0, state, dt, drivingForces);
        end
        
        function [convergence, values, names] = checkConvergence(model, problem, varargin)
            [convergence, values, names] = checkConvergence(model.pressureModel, problem, varargin{:});
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, forces)
            if model.reconstruct
                dispif(model.verbose, 'Pressure converged, reconstructing velocity field...');
                timer = tic();
                [state, report_ms] = model.reconstructFluxes(state0, state, dt, forces);
                dispif(model.verbose, ' Used %s \n', formatTimeRange(toc(timer))); 
            else
                report_ms = [];
            end
            if isfield(state, 'report')
                state = rmfield(state, 'report');
            end
            [state, report] = model.pressureModel.updateAfterConvergence(state0, state, dt, forces);
            report.reconstruction = report_ms;
        end
    
        function [state, varargout] = updateState(model, state, varargin)
            varargout = cell(1, nargout-1);
            p0 = state.pressure;
            [state, varargout{:}] = model.pressureModel.updateState(state, varargin{:});
            if model.coarseScaleConvergence
                part = model.multiscaleSolver.coarsegrid.partition;
                state.dpRel = accumarray(part, abs(state.dpRel))./accumarray(part, 1);
            end
        end
        
        function varargout = getVariableField(model, varargin)
            varargout = cell(1, nargout);
            [varargout{:}] = model.pressureModel.getVariableField(varargin{:});
        end
        
        % Utility functions 
        function [state, report] = reconstructFluxes(model, state0, stateConverged, dt, forces)
            CG = model.multiscaleSolver.coarsegrid;
            nc = CG.parent.cells.num;
            
            timer = tic();
            
            propPressure = stateConverged.pressure;
            
            if isa(model.pressureModel, 'PressureCompositionalModel')
                extra = {'computeFlash', false};
            else
                extra = {};
            end
            problem = model.pressureModel.getEquations(state0, stateConverged, dt, forces,...
                'iteration', inf, 'staticwells', true, extra{:});
            if isa(problem, 'PressureReducedLinearSystem')
                A = problem.A;
                b = problem.b;
            else
                A = problem.equations{1}.jac{1};
                b = -problem.equations{1}.val;
            end
            % Single pass of multiscale solver
            solver = model.multiscaleSolver;
            its = solver.maxIterations;
            solver.maxIterations = 0;

            dpCoarse = solver.solveLinearSystem(A, b);
            
            solver.maxIterations = its;
            
            [dpFlux, t_solve] = reconstructPressure(CG.partition, dpCoarse, A, b);
            
            % Update with reconstruction and prolongation increments
            [stateCoarse, stateFlux] = deal(stateConverged);
            stateCoarse.pressure = propPressure  + dpCoarse(1:nc);
            stateFlux.pressure =  propPressure + dpFlux(1:nc);
            
            % Store fluxes and use it the corect places
            stateFlux   = setFluxes(model, state0, stateFlux,   dt, forces, propPressure);
            stateCoarse = setFluxes(model, state0, stateCoarse, dt, forces, propPressure);
            
            flux = stateFlux.flux;
            flux(CG.faces.fconn, :) = stateCoarse.flux(CG.faces.fconn, :);
            
            state = stateConverged;
            % Take fluxes from wellSol in reconstructed state
            state.wellSol = stateFlux.wellSol;
            state.flux = flux;
            
            % Use property pressure for transport
            state.pressure = propPressure;

            report.reconstructionTime = toc(timer);
            report.reconstructionSolver = t_solve;
        end
        
        function state_flux = setFluxes(model, state0, state, dt, forces, propsPressure)
            [~, state_flux] = model.pressureModel.getEquations(state0, state, dt, forces, ...
                                       'ResOnly', true, 'iteration', inf, ...
                                       'propsPressure', propsPressure, 'staticwells', true);

        end
         
        function state = validateState(model, state)
            state = model.pressureModel.validateState(state);
        end
        
        function [model, state] = updateForChangedControls(model, state, forces)
            [model.pressureModel, state] = model.pressureModel.updateForChangedControls(state, forces);
        end
        
        function [model, state] = prepareReportstep(model, varargin)
            [model.pressureModel, state] = model.pressureModel.prepareReportstep(varargin{:});
        end

        function [model, state] = prepareTimestep(model, varargin)
            [model.pressureModel, state] = model.pressureModel.prepareTimestep(varargin{:});
        end
        function model = validateModel(model, varargin)
            model.pressureModel = model.pressureModel.validateModel(varargin{:});
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
