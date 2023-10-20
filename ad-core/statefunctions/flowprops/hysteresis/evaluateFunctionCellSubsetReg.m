function v = evaluateFunctionCellSubsetReg(prop, fn, regions, varargin)
%
% Based on evaluateFunctionCellSubset within StateFunction. This is used
% here in order to avoid modifying FlowPropertyFunctions (which would need
% to pass different region indicators for different classes).
%
% Used by:
% HystereticRelativePermeability
%
if iscell(fn)
    % We have multiple regions and have to evaluate for each one
    nc = size(regions, 1);
    isCell = cellfun(@(x) numelValue(x) == nc, varargin);
    assert(~isempty(regions))
    [sample, isAD] = getSampleAD(varargin{:});
    v = zeros(numel(regions), 1);
    if isAD
        v = prop.AutoDiffBackend.convertToAD(v, sample);
    end
    for reg = 1:numel(fn)
        if min(regions) > 1 % rock.regions.imbibition passed as regions
            act = regions == reg + (min(regions)-1);
        else
            act = regions == reg;
        end
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