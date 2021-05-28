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

    properties (Access = protected)
        outputRange = [-inf, inf];
        allowInf = false;
        allowNaN = false;
        isSingleRegion = false; % There are multiple regions in fluid, but only a single one is used in practice.
    end
    
    methods
        function prop = StateFunction(model, varargin)
            if nargin > 0
                assert(isa(model, 'PhysicalModel'), 'Model must be derived from PhysicalModel base class');
                prop.AutoDiffBackend = model.AutoDiffBackend;
                if mod(numel(varargin), 2) == 1
                    % Odd number of entries means we specified regions
                    prop.regions = varargin{1};
                    varargin = varargin(2:end);
                end
                prop = merge_options(prop, varargin{:});
            end
            prop.isSingleRegion = isempty(prop.regions) || numel(unique(prop.regions)) == 1;
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
        
        function v = evaluateFluid(prop, model, fluidFunctionName, varargin)
            % Evaluate a named function on the fluid struct on the entire
            % domain with specified function arguments
            fn = model.fluid.(fluidFunctionName);
            v = prop.evaluateFunctionOnDomainWithArguments(fn, varargin{:});
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
            if iscell(fn)
                local_region = prop.regions;
                if isempty(local_region)
                    error(['The function provided is a cell array, which', ...
                           ' indicates that multiple regions are present.', ...
                           ' This instance of %s has empty .regions.', ...
                           ' An region must be provided when', ...
                           ' the input function is a cell array'], class(prop));
                end
                if prop.isSingleRegion
                    lr = local_region(1);
                    v = fn{lr}(varargin{:});
                else
                    if isnumeric(subset) || islogical(subset)
                        local_region = local_region(subset);
                    end
                    % We have multiple regions and have to evaluate for each
                    % subregion
                    nc = size(prop.regions, 1);
                    isCell = cellfun(@(x) numelValue(x) == nc, varargin);
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
            % same StateFunctionGroupings class)
            varargout = cell(nargout, 1);
            s = struct(state.(prop.structName));
            for i = 1:nargout
                v = s.(varargin{i});
                v = expandIfUniform(v);
                varargout{i} = v;
            end
        end
        
        function varargout = getEvaluatedExternals(prop, model, state, varargin)
            % Get evaluated values from external dependencies (from other
            % StateFunctionGroupings class)
            varargout = cell(nargout, 1);
            for i = 1:nargout
                nm = varargin{i};
                isSub = nm == '.';
                if any(isSub)
                    % Syntax GroupingName.StateFunctionName
                    pos = find(isSub);
                    group = nm(1:(pos-1));
                    nm = nm((pos+1):end);
                else
                    pos = arrayfun(@(x) strcmp(x.name, nm), prop.externals);
                    group = prop.externals(pos).grouping;
                end
                if strcmp(group, 'state')
                    v = model.getProp(state, nm);
                else
                    struct_name = model.(group).getStateFunctionContainerName();
                    s = struct(state.(struct_name));
                    v = s.(nm);
                end
                v = expandIfUniform(v);
                varargout{i} = v;
            end
        end
        
        function validateOutput(prop, output, level)
            % Check the outputs at runtime, with more expensive checks
            % occuring for higher levels
            if level == 0
                return
            end
            if isstruct(output)
                % Do nothing. Currently not checking structs
                return;
            end
            % Wrap in cell if not already wrapped
            if ~iscell(output)
                output = {output};
            end
            checkJacobian = level > 2;
            for i = 1:numel(output)
                o = output{i};
                if ~prop.allowNaN
                    prop.validateOutputInternal(o, i, level, 'NaN values', @isnan, checkJacobian);
                end
                if ~prop.allowInf
                    prop.validateOutputInternal(o, i, level, 'Inf values', @isinf, checkJacobian);
                end
                lower = prop.outputRange(1);
                upper = prop.outputRange(2);
                if isfinite(lower)
                    str = sprintf('Value less than %f', lower);
                    prop.validateOutputInternal(o, i, level, str, @(x) x < lower, false);
                end
                if isfinite(upper)
                    str = sprintf('Value larger than %f', upper);
                    prop.validateOutputInternal(o, i, level, str, @(x) x > upper, false);
                end
            end
        end
    end
    
    methods (Access = protected)
        function validateOutputInternal(prop, o, columnIndex, level, descr, test, checkJacobian)
            v = value(o);
            if isnumeric(v)
                bad = test(v);
                numbad = sum(bad);
                if numbad > 0
                    nv = numel(v);
                    warning(['Bad output: %s in %d of %d (%1.2f%%) entries for', ...
                             ' state function ''%s'' in column %d'],...
                            descr, numbad, nv, numbad/nv, class(prop), columnIndex);
                    if level > 2
                        tmp = find(bad);
                        nl = min(numbad, level);
                        fprintf('%d first bad entries:\n', nl);
                        for j = 1:nl
                            index = tmp(j);
                            fprintf('%d: %f\n', index, v(index));
                        end
                    end
                end
                if checkJacobian && isa(o, 'ADI')
                    % Check AD values too.
                    for jacNo = 1:numel(o.jac)
                        bad = test(nonzeros(sparse(o.jac{jacNo})));
                        numbad = sum(bad);
                        if numbad > 0
                            warning(['Bad Jacobian: %s in %d entries for', ...
                                     ' state function ''%s'' column %d Jacobian %d'],...
                                    descr, numbad, class(prop), columnIndex, jacNo);
                        end
                    end
                end
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
