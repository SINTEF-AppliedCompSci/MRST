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
        function problem = LinearizedProblem(equations, types, names, primary, state, dt)
            if nargin < 6
                dt = 0;
            end
            problem.equations = equations;
            problem.types = types;
            problem.equationNames = names;
            problem.primaryVariables = primary;
            problem.A = [];
            problem.b = [];
            problem.dt = dt;
            problem.iterationNo = nan;
            
            problem.state = state;
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
