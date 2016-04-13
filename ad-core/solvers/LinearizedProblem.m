classdef LinearizedProblem
% A linearized problem within a non-linear iteration
properties
    % Cell array of the equations. Can be either doubles, or more typically
    % ADI objects.
    equations
    % Equal length to number of equations, with strings indicating their
    % types (common types: cell for cell variables, well for well equations
    % etc).
    types
    % Equal length to number of equations, giving them unique names for
    % readability.
    equationNames
    % The primary variables used to compute the problem.
    primaryVariables
    % Linear system after assembling Jacobians
    A
    % Right hand side for linearized system
    b
    % The state used to produce the equations
    state
    % The time step lengt
    dt
    % Optionally, the nonlinear iteration number corresponding to this
    % problem.
    iterationNo
    drivingForces
end

methods
    function problem = LinearizedProblem(varargin)
        if nargin == 5 || nargin == 6
            % Call syntax:
            % LinearizedProblem(eqs, eqtypes, eqnames, primaryvars, state, dt)
            eqs         = varargin{1};
            eqtypes     = varargin{2};
            eqNames    = varargin{3};
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
    function [equations, types, names] = checkInputs(problem, equations, types, names)
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
       % Assemble the Jacobian and right hand side and store them if
       % they aren't already created
       if isempty(problem.A)
           % Ignore empty equations
           iseq = cellfun(@(x) ~isempty(x), problem.equations);
           eqs = cat(problem.equations{iseq});
           problem.A = -eqs.jac{1};
           problem.b = eqs.val;
       end
    end
    
    % --------------------------------------------------------------------%
    function problem = clearSystem(problem)
        % Reset A/b
        problem.A = [];
        problem.b = [];
    end
    
    % --------------------------------------------------------------------%
    function [A, b] = getLinearSystem(problem)
        % Get problem suitable for standard linear solvers
        problem = problem.assembleSystem();
        A = problem.A;
        b = problem.b;
    end
    
    % --------------------------------------------------------------------%
    function values = norm(problem, varargin)
        % Overload norm for convergence testing
        values = cellfun(@(x) norm(double(x), varargin{:}), problem.equations);
    end
    
    % --------------------------------------------------------------------%
    function n = numel(problem)
        % Numel gives us the number of equations
        n = numel(problem.equations);
    end
    
    % --------------------------------------------------------------------%
    function problem = reorderEquations(problem, newIndices)
        % Reorder equations to new ordering
        assert(numel(problem) == numel(newIndices));
        
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
        % Get number of variables for one or more equations. Single
        % input argument defaults to all equations.
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
    
    % --------------------------------------------------------------------%
    function [problem, eliminatedEquation] = eliminateVariable(problem, variable)
        % Eliminate a variable from the problem using the equation with
        % the same index.
        if isa(variable, 'char')
            n = find(problem.indexOfEquationName(variable));
        elseif isnumeric(variable)
            n = variable;
            if islogical(n)
                n = find(n);
            end
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
        % Eliminate the non-cell variables first
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
        % Reduced problem (problem resulting from a call to
        % reduceToSingleVariableType)

        % originalProblem is the problem before reduction

        % increments reduced is the cell array of increments from the
        % solution of the reduced problem

        % eliminated is the eliminated equations.

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
