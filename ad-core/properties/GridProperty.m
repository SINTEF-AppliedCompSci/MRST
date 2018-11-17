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

        function v = evaluateFunction(prop, fn, varargin)
            if iscell(fn)
                nc = size(prop.regions, 1);
                isCell = cellfun(@(x) numel(double(x)) == nc, varargin);
                assert(~isempty(prop.regions))
                [sample, isAD] = getSampleAD(varargin{:});
                v = zeros(numel(prop.regions), 1);
                if isAD
                    v = prop.AutoDiffBackend.convertToAD(v, sample);
                end
                for reg = 1:numel(fn)
                    act = prop.regions == reg;
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

        function v = evaluateFunctionSingleRegion(prop, fn, region_index, varargin)
            if iscell(fn)
                v = fn{region_index}(varargin{:});
            else
                v = fn(varargin{:});
            end
        end
    end
end