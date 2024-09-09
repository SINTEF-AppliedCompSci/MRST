function v = evaluateFunctionCellSubsetReg(prop, fn, regions, varargin)
%
% Based on evaluateFunctionCellSubset within StateFunction. This is used
% here in order to avoid modifying FlowPropertyFunctions (which would need
% to pass different region indicators for different classes).
%
% Used by:
% HystereticRelativePermeability

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
