function hash = obj2hash(obj, varargin)
%Compute hash value (md5 checksum) of an object, including classes.
%
% SYNOPSIS:
%   hash = obj2hash(obj)
%   hash = obj2hash(obj, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function computes the hash value (md5 checksum) of an object.
%   Useful for comparing model parameters and uniquely identify
%   simulation setups.
%
% REQUIRED PARAMETERS:
%   obj - object to be hashed
%
% OPTIONAL PARAMETERS:
%   skip    - Field names (if obj is a struct) or property names (if obj is a
%             class) to be skipped in the hash computation.
%             Default value: false
%   maxSize - Maximum number of non-zero elements in obj to be hashed if
%             obj is numeric, logical or a character array. If obj has more
%             nnzs, the function will simply hash the class(obj).
%             Default value: 1e8
%   verbose - Enable output.  Default value dependent upon global verbose
%             settings of function 'mrstVerbose'.
%
% RETURNS:
%   hash - string of 32 characters with the hexadecimal md5 checksum of
%          obj. Specifically, if obj is
%               * numerical, logical or a character array, the hash equals
%                 the md5 checksum of obj (see note on maxSize above).
%               * a cell array, the hash equals the combined hash of the
%                 hash of each element of obj.
%               * a struct, the hash equals the combined hash of the
%                 hash of each of the fields of obj (see note on skipped
%                 field names above).
%               * none of the above, the hash equals the combined hash of
%                 the hash of each of the properties of obj (see note on
%                 skipped field names above). If obj has no public
%                 properties, the function will fall back to hasing
%                 class(obj)
%
% EXAMPLES:
%   obj  = struct('f1', rand(10), 'f2', {{'yes', 'no', struct('a', 2)}});
%   hash = obj2hash(obj);
%
%   G     = computeGeometry(cartGrid([10,10], [100, 100]));
%   rock  = makeRock(G, 100*milli*darcy, 0.1);
%   fluid = initSimpleADIFluid();
%   model = GenericBlackOilModel(G, rock, fluid);
%   hash  = obj2hash(model);
%
% SEE ALSO:
%   `md5sum`

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

    opt = struct('skip'   , {{}}       , ...
                 'maxSize', 1e8        , ...
                 'verbose', mrstVerbose);
    opt = merge_options(opt, varargin{:});
    
    % Avoid empty hash values
    if isempty(obj), obj = '[]'; end
    
    if isnumeric(obj) || islogical(obj) || ischar(obj)
        % Object is a type we can hash directly
        if nnz(obj) > opt.maxSize
            if opt.verbose
                warning(['Object has too many nonzero elements (%d), '  , ...
                         'falling back to hashing class name'], nnz(obj));
            end
            obj = class(obj);
        end
        hash = md5sum(obj);
        return
    elseif iscell(obj)
        % Compute hash of each cell array element and hash combined result
        hash = cellfun(@(o) obj2hash(o, varargin{:}), ...
                       obj, 'UniformOutput', false  );
        hash = combineHashes(hash);
        return
    elseif isstruct(obj)
        % Object is a struct. We will have to loop through all fields
        names = fieldnames(obj);
    else
        % Object is not of a type we can hash directly
        names = properties(obj);
    end
    
    if isempty(names)
        % Object has no public properties. Hash the class name instead
        hash = obj2hash(class(obj));
        return;
    end
    
    % Exclude any fields we have asked to skip
    names = setdiff(names, opt.skip);
    hash  = cell(1, numel(names));
    % Loop through all fields/properties and compute hash
    for i = 1:numel(names)
        % Avoid protected (returned by properties(obj) in Octave)
        if ~isprop(obj, names{i}), continue; end
        prop    = obj.(names{i});
        hash{i} = obj2hash(prop, varargin{:});
    end
    hash = hash(~cellfun(@isempty, hash));
    % Hash the combined result
    hash = combineHashes(hash);
    
end

%-------------------------------------------------------------------------%
function hash = combineHashes(hashes)
    hash = strjoin(hashes, '_');
    hash = md5sum(hash);
end