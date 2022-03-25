function out = compareStructs(in1, in2, varargin)
%Compare fields of two structs
%
% SYNOPSIS:
%   out = compareStructs(in1, in2)
%   out = compareStructs(in1, in2, 'pn1', pv1, ...)
%
% PARAMETERS:
%   in1, in2 - struct to be compared
%
% OPTIONAL ARGUMENTS:
%   'fun'            - Function for comparing fields of the two structs.
%                      For a field 'name', let v1 = in1.name and v2 =
%                      in2.name. Then, the fields will be compared as
%                           out.name = opt.fun(v1) - opt.fun(v2).
%                      Default value is @abs.
%
%   'relative'       - Boolean indicating if the realtive difference should be
%                      computed. If true,
%                           out.name = (opt.fun(v1) - opt.fun(v2))./ ...
%                               max(max(opt.fun(v1), opt.fun(v2)), 1e-10);
%                      Default value is false.
%
%   'skip'           - Cell array of field names to be skipped in the
%                      comparison. Default value is empty.
%
%   'includeStructs' - Boolean indicating if fields that are themselves
%                      structs should be included in the comparison.
%
%   'verbose'        - Enable output. Default value dependent upon global
%                      verbose settings of function 'mrstVerbose'.
%
%
% RETURNS:
%   out - struct with one field for each of the fields in in1 and in2 that
%         has been compared.
%
% EXAMPLE:
%   in1 = struct('a', 1:10);
%   in2 = struct('a', 2:11);
%   out = compareStructs(in1, in2, 'relative', true);
%

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('fun'           , @abs         , ...
                 'relative'      , false        , ...
                 'skip'          , {{}}         , ...
                 'includeStructs', true         , ...
                 'verbose'       , mrstVerbose());
    opt = merge_options(opt, varargin{:});
    
    if numel(in1) > 1
        n1 = numel(in1);
        n2 = numel(in2);
        out = [];
        for i = 1:min(n1, n2)
            o   = compareStructs(in1(i), in2(i), varargin{:});
            out = [out; o]; %#ok
            return;
        end
    end
    
    assert(isstruct(in1) && isstruct(in2), 'Both inputs must be structs');
    
    keep1 = getValidFields(in1, opt);
    keep2 = getValidFields(in2, opt);
    
    names1 = fieldnames(in1); names1 = names1(keep1);
    names2 = fieldnames(in2); names2 = names2(keep2);
    
    names   = intersect(names1, names2);
    names12 = union(names1, names2);
    
    unhandled = setdiff(names12, names);
    if ~isempty(unhandled) && opt.verbose
        warning('Found %d fields not common to both structs', numel(unhandled));
    end
    
    out = cell2struct(num2cell(nan(numel(names12),1)), names12);
    
    for i = 1:numel(names)
        name = names{i};
        v1 = in1.(name);
        v2 = in2.(name);
        if isstruct(v1) && opt.includeStructs()
            assert(isstruct(v2))
            out.(name) = compareStructs(v1, v2, varargin{:});
        else
            if iscell(v1)
                assert(iscell(v2))
                v = cellfun(@(v1,v2) compareValues(v1, v2, opt), ...
                                           v1, v2, 'UniformOutput', false);
            else
                v = compareValues(v1, v2, opt);
            end
            out.(name) = v;
        end
    end
    
end

%-------------------------------------------------------------------------%
function keep = getValidFields(in, opt)
% Get valid candidate fields for comparison in struct

    v = struct2cell(in);
    n = numel(v);
    keep = false(n, 1);
    keepfun = @(v) isnumeric(v) || islogical(v);
    if opt.includeStructs
        keepfun = @(v) keepfun(v) || isstruct(v);
    end
    for i = 1:n
        if iscell(v{i})
            keep(i) = all(cellfun(@(v) keepfun(v), v{i}));
        else
            keep(i) = keepfun(v{i});
        end 
    end
    if ~isempty(opt.skip)
        keep = keep & cellfun(@(n) ~any(strcmpi(n, opt.skip)), fieldnames(in));
    end
    
end

%-------------------------------------------------------------------------%
function v = compareValues(v1, v2, opt)
% Compute difference of opt.fun(v1) and opt.fun(v2), where opt.fun
% typically is a norm

    v = opt.fun(v1 - v2);
    if opt.relative
        v = v./max(max(opt.fun(v1), opt.fun(v2)), 1e-10);
    end
    
end
