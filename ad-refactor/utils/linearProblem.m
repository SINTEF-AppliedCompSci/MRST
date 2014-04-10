classdef linearProblem
    % A linear subproblem within a non-linear iteration
    properties
        equations
        types
        equationNames
        primaryVariables
        A
        b
    end
    
    methods
        function problem = linearProblem(equations, types, names, primary)
            problem.equations = equations;
            problem.types = types;
            problem.equationNames = names;
            problem.primaryVariables = primary;
            problem.A = [];
            problem.b = [];
        end
        
        function problem = assembleSystem(problem)
           % Assemble the Jacobian and right hand side and store them if
           % they aren't already created
           if isempty(problem.A)
               eqs = cat(problem.equations{:});
               problem.A = -eqs.jac{1};
               problem.b = eqs.val;
           end
        end
        
        function values = norm(problem, varargin)
            % Overload norm for convergence testing
            values = cellfun(@(x) norm(x.val, varargin{:}), problem.equations);
        end
    end
end