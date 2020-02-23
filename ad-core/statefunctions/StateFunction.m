classdef StateFunction
    % Class for gridded domain properties
    properties
        regions
        AutoDiffBackend
        dependencies = {};
        externals = [];
        structName
        label
    end

    methods
        function prop = StateFunction(model, regions, varargin)
            if nargin > 0
                assert(isa(model, 'PhysicalModel'), 'Model must be derived from PhysicalModel base class');
                prop.AutoDiffBackend = model.AutoDiffBackend;
                if nargin > 1
                    prop.regions = regions;
                end
            end
        end

        function value = evaluateOnDomain(prop, model, state)
            % Given state, evaluate the canonical representation for the
            % current model.
            error('Base class should not be evaluated')
        end

        function prop = dependsOn(prop, varargin)
            % Document dependencies and external dependencies
            prop = addPropertyDependence(prop, varargin{:});
        end
        
        function prop = resetDependencies(prop, resetInternal, resetExternal)
            if nargin < 3
                resetExternal = true;
            end
            if nargin < 2
                resetInternal = true;
            end
            if resetInternal
                prop.dependencies = {};
            end
            if resetExternal
                prop.externals = [];
            end
        end
        
        function v = evaluateFunctionOnDomainWithArguments(prop, fn, varargin)
            % Evaluate function handle on the entire domain, with specified
            % input arguments
            v = prop.evaluateFunctionCellSubset(fn, ':', varargin{:});
        end

        function v = evaluateFunctionSingleRegionWithArguments(prop, fn, region_index, varargin)
            % Evaluate function on a single specific region, with specified
            % input arguments
            assert(region_index <= numel(fn), 'Region index exceeds maximum number of regions.');
            if iscell(fn)
                v = fn{region_index}(varargin{:});
            else
                v = fn(varargin{:});
            end
        end
        
        function v = evaluateFunctionCellSubset(prop, fn, subset, varargin)
            % Evaluate specific function on a given subset
            if ischar(subset) && strcmp(subset, ':')
                local_region = prop.regions;
            else
                local_region = prop.regions(subset);
            end
            if iscell(fn)
                % We have multiple regions and have to evaluate for each
                % subregion
                nc = size(prop.regions, 1);
                isCell = cellfun(@(x) numelValue(x) == nc, varargin);
                assert(~isempty(prop.regions))
                [sample, isAD] = getSampleAD(varargin{:});
                v = zeros(numel(local_region), 1);
                if isAD
                    v = prop.AutoDiffBackend.convertToAD(v, sample);
                end
                for reg = 1:numel(fn)
                    act = local_region == reg;
                    arg = varargin;
                    carg = cellfun(@(x) x(act), arg(isCell), 'UniformOutput', false);
                    [arg{isCell}] = carg{:};
                    if any(act)
                        v(act) = fn{reg}(arg{:});
                    end
                end
            else
                v = fn(varargin{:});
            end
        end
        
        function property = subset(property, cell_subset)
            % Take a subset-property (e.g. for extracting the function on a
            % local domain)
            if ~isempty(property.regions)
                property.regions = property.regions(cell_subset);
            end
        end
        
        function varargout = getEvaluatedDependencies(prop, state, varargin)
            % Get evaluated values from local dependencies (belonging to
            % same PropertyFunctions grouping)
            varargout = cell(nargout, 1);
            s = struct(state.(prop.structName));
            for i = 1:nargout
                v = s.(varargin{i});
                v = expandIfUniform(v);
                varargout{i} = v;
            end
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
