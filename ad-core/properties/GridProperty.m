classdef GridProperty
    properties
        regions
        AutoDiffBackend
    end

    methods
        function prop = GridProperty(backend, regions, varargin)
            prop.AutoDiffBackend = backend;
            prop.regions = regions;
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
            local_region = prop.regions(subset);
            if iscell(fn)
                % We have multiple regions and have to evaluate for each
                % subregion
                nc = size(prop.regions, 1);
                isCell = cellfun(@(x) numel(double(x)) == nc, varargin);
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
    end
end