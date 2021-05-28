classdef StateFunctionDependent
    % A virtual class used for multiple inheritance. If a class inherits
    % from this class, it can use the function function dependsOn to
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
