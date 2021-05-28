classdef LinearizedProblem
    % A linearized problem within a non-linear iteration
    %
    % SYNOPSIS:
    %   problem = LinearizedProblem(primvars, state)
    %   problem = LinearizedProblem(primvars, state, dt)
    %   problem = LinearizedProblem(eqs, eqtypes, eqnames, primvars, state)
    %   problem = LinearizedProblem(eqs, eqtypes, eqnames, primvars, state, dt)
    %
    % DESCRIPTION:
    %   A class that implements storage of an instance of a linearized problem
    %   discretized using AD. This class contains the residual equations
    %   evaluated for a single state along with meta-information about the
    %   equations and the primary variables they are differentiated with
    %   respect to.
    %
    %   A linearized problem can be transformed into a linear system and solved
    %   using `LinearSolverAD`-derived subclasses, given that the number of
    %   equations matches the number of primary variables.
    %
    %   In particular, the class contains member functions for:
    %
    %     - assembling a linear system from the Jacobian block matrices stored
    %       for each individual (continuous) equation
    %     - appending/prepending additional equations
    %     - using a block-Gaussian method to eliminate individual variables or
    %       all variables that are not of a specified type
    %     - recovering increments corresponding to variables that have
    %       previously been eliminated
    %     - computing the norm of each residual equation
    %
    %   as well as a number of utility functions for sanity checks, quering of
    %   indices and the number of equations, clearing the linear system, etc.
    %
    % SEE ALSO:
    %   `LinearSolverAD`, `PhysicalModel`

properties
    equations % Cell array of the equations. Can be doubles, or more typically ADI objects.
    types % Cell array of equal length to number of equations, with strings indicating their types 
    % (common types: cell for cell variables, well
    % for well equations etc). Note that these types are available to any
    % linear solver, which may use them to construct appropriate solver
    % strategies or preconditioners.
    equationNames% Cell array of equal length to number of equations, giving them unique names
    % that are used e.g. when printing when convergence reports or solving
    % the linear system.
    primaryVariables % Cell array containing the names of the primary variables.
   
    A % Linear system after assembling Jacobians. 
    % Will be empty until the member function assembleSystem has been called.
    b % Right hand side for linearized system. See A.
    state % The problem state used to produce the equations
    dt % The time step length (Not relevant for all problems)
    iterationNo % Nonlinear iteration number corresponding to problem.
    drivingForces % The driving forces struct used to generate the equations.
end

methods
    function problem = LinearizedProblem(varargin)
        if nargin == 5 || nargin == 6
            % Call syntax:
            % LinearizedProblem(eqs, eqtypes, eqnames, primaryvars, state, dt)
            eqs         = varargin{1};
            eqtypes     = varargin{2};
            eqNames     = varargin{3};
            primaryVars = varargin{4};
            modelState  = varargin{5};

        if nargin == 5
            timestep = 0;
        else
            timestep = varargin{6};
        end

        elseif nargin == 2 || nargin == 3
            % Call syntax:
            % LinearizedProblem(primaryVariables, state, dt)
            [eqs, eqtypes, eqNames] = deal({});

            primaryVars = varargin{1};
            modelState = varargin{2};
            if nargin == 2
                timestep = 0;
            else
                timestep = varargin{3};
            end
        else
            error('Bad number of input arguments');
        end
        [eqs, eqtypes, eqNames] = problem.checkInputs(eqs, eqtypes, eqNames);
        if ischar(primaryVars)
            primaryVars = {primaryVars};
        end

        problem.equations = eqs;
        problem.types = eqtypes;
        problem.equationNames = eqNames;
        problem.primaryVariables = primaryVars;

        assert(all([numel(eqtypes), numel(eqNames)] == numel(eqs)), ...
            'Inconsistent number of types/names/equation numbers.');

        problem.A = [];
        problem.b = [];
        problem.dt = timestep;
        problem.iterationNo = nan;

        problem.state = modelState;
    end
    
    % --------------------------------------------------------------------%
    function problem = prependEquations(problem, equations, types, names)
        % Add one or more equations to the end of the current list of
        % equations.
        [equations, types, names] = problem.checkInputs(equations, types, names);

        problem.equations = [equations, problem.equations];
        problem.types     = [types, problem.types];
        problem.equationNames     = [names, problem.equationNames];

        % Reset linear system
        problem.A = [];
        problem.b = [];
    end
    
    % --------------------------------------------------------------------%
    function problem = appendEquations(problem, equations, types, names)
        % Add one or more equations to the beginning of the current
        % list of equations.
        [equations, types, names] = checkInputs(problem, equations, types, names);
        problem.equations = [problem.equations, equations];
        problem.types     = [problem.types, types];
        problem.equationNames     = [problem.equationNames, names];
        % Reset linear system
        problem.A = [];
        problem.b = [];
    end
    
    % --------------------------------------------------------------------%
    function [equations, types, names] = checkInputs(problem, equations, types, names) %#ok
        if iscell(equations)
            % For multiple inputs, we want them to be of the same
            % length
            assert(all([numel(types), numel(names)] == numel(equations)), ...
                'Inconsistent number of types/names/equation numbers.');
        else
            equations = {equations};
            types = {types};
            names = {names};
        end
    end
    
    % --------------------------------------------------------------------%
    function problem = assembleSystem(problem)
       % Assemble the linear system from the individual Jacobians and
       % residual functions. Note the return parameter, as problem is not a
       % handle class.
       if isempty(problem.A)
           % Ignore empty equations
           iseq = cellfun(@(x) ~isempty(x), problem.equations);
           eqs = combineEquations(problem.equations{iseq});
           if isa(eqs, 'ADI')
               problem.A = eqs.jac{1};
               problem.b = -eqs.val;
           else
               problem.b = -eqs;
           end
       end
    end
    % --------------------------------------------------------------------%
    function [result, report] = processResultAfterSolve(problem, result, report) %#ok
        % Do nothing
    end
    
    % --------------------------------------------------------------------%
    function problem = clearSystem(problem)
        % Remove pre-assembled linear system. Typically used if the
        % equations are manually changed, and one does not want to evaluate
        % inconsistent linear systems.
        problem.A = [];
        problem.b = [];
    end
    
    % --------------------------------------------------------------------%
    function [A, b] = getLinearSystem(problem)
        % Get the linear system from the individual Jacobians and residual
        % equations in a format suitable for general linear solvers. If the
        % equations are not already assembled, this will result in a call
        % to assembleSystem. If you will retrieve the linear system
        % multiple times, it is better to call assembleSystem beforehand.
        problem = problem.assembleSystem();
        A = problem.A;
        b = problem.b;
    end
    
    % --------------------------------------------------------------------%
    function values = norm(problem, varargin)
        % Get the norm of each equation. This is an overloaded function.
        values = cellfun(@(x) norm(value(x), varargin{:}), problem.equations);
    end
    
    % --------------------------------------------------------------------%
    function n = numeq(problem)
        % Get the number of distinct equations. Note that an equation here
        % can refer to multiple subequations of the same type (e.g. a 100
        % cell problem where the equations are oil and water conservation
        % laws will have two equations in total, with each having 100
        % sub-entries.
        n = numel(problem.equations);
    end
    
    % --------------------------------------------------------------------%
    function problem = reorderEquations(problem, newIndices)
        % Reorder equations based on a set of indices.
        % ARGUMENTS:
        % - newIndices: numeq(problem) long array of new indices into the
        %               equations.
        % OUTPUT:
        % - problem   : Problem with re-numbered equations.
        assert(numeq(problem) == numel(newIndices));
        
        problem.equations = problem.equations(newIndices);
        if ~isempty(problem.types)
            problem.types = problem.types(newIndices);
        end
        if ~isempty(problem.equationNames)
            problem.equationNames = problem.equationNames(newIndices);
        end
    end
    
    % --------------------------------------------------------------------%
    function varnum = getEquationVarNum(problem, n)
        % Get number of subequations for one or more equations. A single
        % equation is a collection of a number of equations grouped by the
        % constructor. Typically, all subequations should be of the same
        % ARGUMENTS:
        %   n -   (OPTIONAL) Indices of the equations for which the number
        %         of subequations is desired. If omitted, the function
        %         returns the number for all equations.
        %
        % RETURNS:
        %   varnum - Array with the number of subequations for the requested
        %            equations.
        if nargin == 1
            n = ':';
        end
        varnum = cellfun(@(x) numel(x.val), problem.equations(n));
    end
    
    % --------------------------------------------------------------------%
    function index = indexOfType(problem, name)
        % Get the index(es) of a type of variable by name.
        index = cellfun(@(x) strcmpi(x, name), problem.types);
    end
    
    % --------------------------------------------------------------------%
    function no = countOfType(problem, name)
        % Count the number of equations that are of a specific type.
        no = sum(cellfun(@(x) strcmpi(x, name), problem.types));
    end
    
    % --------------------------------------------------------------------%
    function index = indexOfPrimaryVariable(problem, name)
        % Get the index of a primary variable by name.
        index = cellfun(@(x) strcmpi(x, name), problem.primaryVariables);
    end
    
    % --------------------------------------------------------------------%
    function index = indexOfEquationName(problem, name)
        % Get the index into the list of equations for a specific name
        index = cellfun(@(x) strcmpi(x, name), problem.equationNames);
    end
    
    function [dvar, dx, problem] = popIncrement(problem, dx, name)
        % Get named/logic position update from a set of increments.
        % Optionally, return dx and problem where the increment has then
        % been removed.
        if ischar(name)
            pos = strcmpi(problem.primaryVariables, name);
        else
            pos = name;
            assert(islogical(pos));
        end
        if any(pos)
            dvar = dx{pos};
            if nargout > 1
                dx = dx(~pos);
                if nargout > 2
                    problem.primaryVariables = problem.primaryVariables(~pos);
                end
            end
        else
            dvar = [];
        end
    end
    
    % --------------------------------------------------------------------%
    function [problem, eliminatedEquation] = eliminateVariable(problem, variable)
        % Eliminate a variable from the problem using the equation with
        % the same index.
        %
        % ARGUMENTS: 
        %   variable - The name of the variable (corresponding to an entry in
        %              problem.names) that is to be eliminated.
        %
        % RETURNS:
        %   problem - Modified problem.
        %   eliminatedEquation - The equation that was eliminated.
        %
        % NOTE:
        %   For non-diagonal matrices, the cost of the equation elimination
        %   can be large.
        if isa(variable, 'char')
            n = find(problem.indexOfEquationName(variable));
        elseif isnumeric(variable)
            n = variable;
        elseif islogical(variable)
            n = find(variable);
        end
        eqs = problem.equations;
        solveInx = setdiff(1:numel(eqs), n);
        eliminatedEquation      = eqs{n};
        for eqNum = solveInx
            for jacNum = solveInx
                if numel(eqs{eqNum}.jac{jacNum}) ~= 0 && numel(eliminatedEquation.jac{jacNum}) ~= 0
                    eqs{eqNum}.jac{jacNum} = eqs{eqNum}.jac{jacNum} - eqs{eqNum}.jac{n}*(eliminatedEquation.jac{n}\eliminatedEquation.jac{jacNum});
                end
            end
            if ~isempty(eqs{eqNum}.val) && ~isempty(eliminatedEquation.val) 
                eqs{eqNum}.val = eqs{eqNum}.val - eqs{eqNum}.jac{n}*(eliminatedEquation.jac{n}\eliminatedEquation.val);
            end
        end

        eqs  = eqs(solveInx);
        for eqNum = 1:numel(eqs)
            eqs{eqNum}.jac = eqs{eqNum}.jac(solveInx);
        end
        problem.equations = eqs;
        problem.equationNames  = problem.equationNames(solveInx);
        problem.types          = problem.types(solveInx);

        % We are implicitly eliminating the variable of the same type
        problem.primaryVariables = problem.primaryVariables(solveInx);
    end
    
    % --------------------------------------------------------------------%
    function [problem, eliminated] = reduceToSingleVariableType(problem, type)
        % Eliminate all equations that are not of a given type
        %
        % ARGUMENTS: 
        %   type - String matching one of problem.types. Equations that DO NOT
        %          have that type will be eliminated.
        %
        % RETURNS:
        %   problem    - Modified problem.
        %   eliminated - Eliminated equations. Can be used to recover the
        %                eliminated values later.
        %
        % NOTE:
        %   Can have a severe cost depending on the sparsity patterns
        %   involved in the different equations. Typical usage is to
        %   eliminate non-cell equations (wells, control equations).
        isCurrent = problem.indexOfType(type);

        % Eliminate all equations that are not of that type
        problem = problem.clearSystem();

        notCellIndex = find(~isCurrent);
        eliminated = cell(numel(notCellIndex), 1);
        elimNames = problem.equationNames(notCellIndex);
        % Gradually peel off problems
        for i = 1:numel(notCellIndex)
            [problem, eliminated{i}] = problem.eliminateVariable(elimNames{i});
        end

    end
    
    % --------------------------------------------------------------------%
    function dx = recoverFromSingleVariableType(reducedProblem, originalProblem, incrementsReduced, eliminated)
        % Recover the increments in primary variables corresponding to
        % equations that have previously been eliminated by for instance
        % "reduceToSingleVariableType".
        % ARGUMENTS:
        %   originalProblem   - The problem before elimination.
        %   incrementsReduced - Increments from the solution of the reduced
        %                       problem.
        %   eliminated        - The eliminated equations.
        %
        % RETURNS:
        %   dx - All increments, as if the `originalProblem` was solved
        %        directly.

        nP = numel(originalProblem);

        type = reducedProblem.types{1};
        current = originalProblem.indexOfType(type);
        notCellIndex = find(~current);
        cellIndex = find(current);
        cellEqNo = numel(cellIndex);

        % Set up storage for all variables, including those we
        % eliminated previously
        dx = cell(nP, 1);

        % Recover non-cell variables
        recovered = false(nP, 1);
        recovered(cellIndex) = true;

        % Put the recovered variables into place
        dx(recovered) = incrementsReduced;

        assert(all(diff(cellIndex) == 1), 'This solver currently assumes that the cell variables comes first!')
        for i = numel(eliminated):-1:1
            pos = notCellIndex(i);
            dVal = recoverVars(eliminated{i}, cellEqNo + 1, dx(recovered));
            dx{pos} = dVal;
            recovered(pos) = true;
        end
    end
end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

