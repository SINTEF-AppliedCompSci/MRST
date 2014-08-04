classdef PhysicalModel
%Base class for physical models
%
% SYNOPSIS:
%   model = PhysicalModel(G, rock, fluid)
%
% DESCRIPTION:
%   Base class for implementing physical models for use with automatic
%   differentiation. This class cannot be used directly.
%
%   A physical model instance contains the functions for getting residuals
%   and jacobians, making a single nonlinear step and verifying
%   convergence. It also contains the functions for updating the state
%   based on the increments found by the linear solver so that the values
%   are physically correct.
%
% REQUIRED PARAMETERS:
%   G     - Simulation grid.
%
%   rock  - Valid rock used for the model.
%
%   fluid - Fluid model used for the model.
%
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   See class properties.
%
% RETURNS:
%   Class instance.
%
% SEE ALSO:
%   ThreePhaseBlackOilModel, TwoPhaseOilWaterModel

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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

    properties
        % Unique identifier for the model
        name
        % The fluid model
        fluid
        % Operators used for construction of systems
        operators
        % Inf norm tolerance for nonlinear iterations
        nonlinearTolerance
        % Grid
        G
        
        % Maximum relative pressure change
        dpMax
        % Maximum relative saturation change
        dsMax
        
        % Water phase present
        water
        % Gas phase present
        gas
        % Oil phase present
        oil
        
        % Names of each component, corresponding to their order in state.s
        componentNames
        
        % Input data used to instantiate the model
        inputdata
    end
    
    methods
        function model = PhysicalModel(G, rock, fluid, varargin) %#ok
            model.dpMax = inf;
            model.dsMax = .2;
            model.nonlinearTolerance = 1e-6;
            model.inputdata = [];
            
            % Physical model
            model.G = G;
            model.fluid = fluid;
            model.componentNames = {'sw', 'so', 'sg'};
            
            model = merge_options(model, varargin{:});
            
            % Base class does not support any phases
            model.water = false;
            model.gas = false;
            model.oil = false;
            
        end
        
        function model = setupOperators(model, G, rock, varargin)
            % Set up divergence/gradient/transmissibility operators
            model.operators = setupSimComp(G, rock, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, drivingForces, dt, varargin) %#ok
            % Get the equations governing the system
            error('Base class not meant for direct use')
        end
        
        function [state, report] = updateState(model, state, dx, drivingForces) %#ok
            % Update state based on non-linear increment
            error('Base class not meant for direct use')
        end
        
        function state = updateAfterConvergence(model, state0, state, drivingForces) %#ok
            % Update state based on non-linear increment after timestep has
            % converged. Defaults to doing nothing since not all models
            % require this.
        end
        
        function [convergence, values] = checkConvergence(model, problem, n)
            % Check and report convergence based on residual tolerances
            if nargin == 2
                n = inf;
            end
            
            values = norm(problem, n);
            convergence = all(values < model.nonlinearTolerance);
            
            if mrstVerbose()
                for i = 1:numel(values)
                    fprintf('%s (%s): %2.2e\t', problem.equationNames{i}, problem.types{i}, values(i));
                end
                fprintf('\n')
            end
        end
        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolve, nonlinsolve, varargin)
            % Make a single linearized timestep
            [problem, state] = model.getEquations(state0, state, dt, drivingForces, varargin{:});
            [convergence, values] = model.checkConvergence(problem);
            
            if ~convergence
                % Get increments for Newton solver
                [dx, ~, linearReport] = linsolve.solveLinearProblem(problem, model);
                
                % Let the non-linear solver decide what to do with the
                % increments to get the best convergence
                dx = nonlinsolve.stabilizeNewtonIncrements(problem, dx);
                
                % Finally update the state. The physical model knows which
                % properties are actually physical.
                [state, updateReport] = model.updateState(state, problem, dx, drivingForces);
                failure = any(cellfun(@(d) ~all(isfinite(d)), dx));
            else
                failure = false;
                [linearReport, updateReport] = deal(struct());
            end
            report = struct('LinearSolver', linearReport, ...
                            'UpdateState',  updateReport, ...
                            'Failure',      failure, ...
                            'Converged',    convergence, ...
                            'Residuals',    values);
        end
        
        function [gradient, result, report] = solveAdjoint(model, solver, getState,...
                                    getObjective, schedule, gradient, itNo, scaling)
            % Solve adjoints
            if nargin == 7
               scaling = struct('rate', 1, 'pressure', 1);
            end
            
            dt_steps = schedule.step.val;
            
            current = getState(itNo);
            before    = getState(itNo - 1);
            dt = dt_steps(itNo);
            
            lookupCtrl = @(step) schedule.control(schedule.step.control(step));
            [~, forces] = model.getDrivingForces(lookupCtrl(itNo));
            problem = model.getEquations(before, current, dt, forces, 'iteration', inf, 'scaling', scaling);
            
            if itNo < numel(dt_steps)
                after    = getState(itNo + 1);
                dt_next = dt_steps(itNo + 1);
                
                [~, forces_p] = model.getDrivingForces(lookupCtrl(itNo + 1));
                problem_p = model.getEquations(current, after, dt_next, forces_p,...
                                    'iteration', inf, 'reverseMode', true, 'scaling', scaling);
            else
                problem_p = [];
            end
            [gradient, result, rep] = solver.solveAdjointProblem(problem_p,...
                                        problem, gradient, getObjective(itNo), model);
            report = struct();
            report.Types = problem.types;
            report.LinearSolverReport = rep;
        end
        
        function [vararg, driving] = getDrivingForces(model, control) %#ok
            % Setup and pass on driving forces
            vararg = {};
            driving = struct('Wells', [], 'bc', [], 'src', []);
            
            if isfield(control, 'W') && ~isempty(control.W)
                vararg = [vararg, 'Wells', control.W];
                driving.Wells = control.W;
            end

            if isfield(control, 'bc') && ~isempty(control.bc)
                vararg = [vararg, 'bc', control.bc];
                driving.bc = control.bc;
            end
            
            if isfield(control, 'src') && ~isempty(control.src)
                vararg = [vararg, 'src', control.src];
                driving.src = control.src;
            end
        end
        
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            index = [];
            switch(lower(name))
                case {'t', 'temperature'}
                    fn = 'T';
                case {'sw', 'water'}
                    index = find(strcmpi(model.componentNames, 'sw'));
                    fn = 's';
                case {'so', 'oil'}
                    index = find(strcmpi(model.componentNames, 'so'));
                    fn = 's';
                case {'sg', 'gas'}
                    index = find(strcmpi(model.componentNames, 'sg'));
                    fn = 's';
                case {'s', 'sat', 'saturation'}
                    index = 1:numel(model.componentNames);
                    fn = 's';
                case {'pressure', 'p'}
                    index = 1;
                    fn = 'pressure';
            end
            
            if isempty(index)
                error('PhysicalModel:UnknownPhase', ...
                    ['Phase ''', name, ''' is not known to this model']);
            end
        end
        
        function p = getProp(model, state, name)
            % Get a property based on the name
            [fn, index] = model.getVariableField(name);
            p = state.(fn)(:, index);
        end
        
        function state = incrementProp(model, state, name, increment)
            % Increment property based on name
            [fn, index] = model.getVariableField(name);
            p = state.(fn)(:, index)  + increment;
            state.(fn)(:, index) = p;
        end
        
        function state = setProp(model, state, name, value)
            % Set property to given value based on name
            [fn, index] = model.getVariableField(name);
            state.(fn)(:, index) = value;
        end
        
        function dv = getIncrement(model, dx, problem, name)
            % Find increment in linearized problem with given name, or
            % output zero if not found
            isVar = problem.indexOfPrimaryVariable(name);
            if any(isVar)
                dv = dx{isVar};
            else
                dv = 0;
            end
        end
        
        function [state, val, val0] = updateStateFromIncrement(model, state, dx, problem, name, relchangemax, abschangemax)
            % Update a state, with optionally a maximum relative change
            % applied.
            if nargin < 6
                relchangemax = inf;
            end
            
            if nargin < 7
                abschangemax = inf;
            end
            
            val0 = model.getProp(state, name);
            if iscell(dx)
                dv = model.getIncrement(dx, problem, name);
            else
                % Numerical value, increment directly and do not safety
                % check that this is a part of the model
                dv = dx;
            end
            biggestChange = max(abs(dv./val0), [], 2);
            relativeChange = min(relchangemax./biggestChange, 1);
            
            biggestChange = max(abs(dv), [], 2);
            absChange = min(abschangemax./biggestChange, 1);
            
            change = min(absChange, relativeChange);
            
            val     = val0 + dv.*repmat(change, 1, size(dv, 2));
            state = model.setProp(state, name, val);

        end
            
        function [isActive, phInd] = getActivePhases(model)
            isActive = [model.water, model.oil, model.gas];
            if nargout > 1
                phInd = find(isActive);
            end
        end
        
        function phNames = getPhaseNames(model)
            tmp = 'WOG';
            active = model.getActivePhases();
            phNames = tmp(active);
        end
    end

    methods (Static)

    end

end

