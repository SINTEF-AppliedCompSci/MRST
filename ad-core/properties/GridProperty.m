classdef GridProperty
    % Class for gridded domain properties
    properties
        regions
        AutoDiffBackend
    end

    methods
        function prop = GridProperty(backend, regions, varargin)
            if nargin > 0
                prop.AutoDiffBackend = backend;
                if nargin > 1
                    prop.regions = regions;
                end
            end
        end

        function value = evaluateOnDomain(prop, model, state)

        end

        function v = evaluateFunctionOnGrid(prop, fn, varargin)
            v = prop.evaluateFunctionCellSubset(fn, ':', varargin{:});
        end

        function v = evaluateFunctionSingleRegion(prop, fn, region_index, varargin)
            assert(region_index <= numel(fn), 'Region index exceeds maximum number of regions.');
            if iscell(fn)
                v = fn{region_index}(varargin{:});
            else
                v = fn(varargin{:});
            end
        end
        
        function v = evaluateFunctionCellSubset(prop, fn, subset, varargin)
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
            if ~isempty(property.regions)
                property.regions = property.regions(cell_subset);
            end
        end
    end
end