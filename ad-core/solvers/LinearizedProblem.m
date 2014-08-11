classdef LinearizedProblem
    % A linear subproblem within a non-linear iteration
    properties
        equations
        types
        equationNames
        primaryVariables
        A
        b
        state
        dt
        iterationNo
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
        
        function problem = clearSystem(problem)
            problem.A = [];
            problem.b = [];
        end
        
        function [A, b] = getLinearSystem(problem)
            problem = problem.assembleSystem();
            A = problem.A;
            b = problem.b;
        end
        
        %%%%%  Overloaded functions  %%%%%
        function values = norm(problem, varargin)
            % Overload norm for convergence testing
            values = cellfun(@(x) norm(double(x), varargin{:}), problem.equations);
        end
        
        function n = numel(problem)
            n = numel(problem.equations);
        end
        %%%%% Utilities %%%%
        function varnum = getEquationVarNum(problem, n)
            if nargin == 1
                n = ':';
            end
            varnum = cellfun(@(x) numel(x.val), problem.equations(n));
        end
        
        function index = indexOfType(problem, name)
            index = cellfun(@(x) strcmpi(x, name), problem.types);
        end
        
        function no = countOfType(problem, name)
            no = sum(cellfun(@(x) strcmpi(x, name), problem.types));
        end
        
        function index = indexOfPrimaryVariable(problem, name)
            index = cellfun(@(x) strcmpi(x, name), problem.primaryVariables);
        end
        
        function index = indexOfEquationName(problem, name)
            index = cellfun(@(x) strcmpi(x, name), problem.equationNames);
        end
        
        %%%%% Linear magic... %%%%%
        function [problem, eliminatedEquation] = eliminateVariable(problem, variable)
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
        
        
    end
end
