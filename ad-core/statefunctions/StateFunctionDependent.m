classdef StateFunctionDependent
    % A virtual class used for multiple inheritance. If a class inherits
    % from this class, it can use the function functionDependsOn to
    % document which member functions depend on specific state functions.
    properties (Access = protected)
        functionDependencies = struct();
    end
    
    methods
        function sfd = StateFunctionDependent()
            % Base constructor
        end
        
        function prop = functionDependsOn(prop, fn, varargin)
            % Document the dependencies of a specific function, either
            % internal or external.
            if ~isfield(prop.functionDependencies, fn)
                prop = prop.clearFunctionDependencies(fn);
            end
            f = prop.functionDependencies;
            f.(fn) = addPropertyDependence(f.(fn), varargin{:});
            prop.functionDependencies = f;
        end
        
        function dep = getFunctionDependencies(prop, fn)
            % Clear dependencies for a specific function
            fd = prop.functionDependencies;
            if isfield(fd, fn)
                dep = fd.(fn);
            else
                dep = struct('dependencies', [], 'externals', []);
            end
        end

        function prop = clearFunctionDependencies(prop, fn)
            % Clear dependencies for a specific function
            prop.functionDependencies.(fn) = struct('dependencies', [], 'externals', []);
        end
    end
end