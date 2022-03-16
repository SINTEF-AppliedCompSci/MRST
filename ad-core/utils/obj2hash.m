function [hash, hashStruct] = obj2hash(obj, varargin)
%Compute hash value (md5 checksum) of an object, including classes.
%
% SYNOPSIS:
%   [hash, hashStruct] = obj2hash(obj)
%   [hash, hashStruct] = obj2hash(obj, 'pn1', pv1, ...)
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
%   skip    - Field names (if obj is a struct) or property names (if obj is
%             a class) to be skipped in the hash computation.
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
%               1) numerical, logical or a character array, the hash equals
%                  the md5 checksum of obj (see note on maxSize above).
%               2) a cell array, the hash equals the combined hash of the
%                  hash of each element of obj.
%               3) a struct or table, the hash equals the combined hash of
%                  the hash of each of the fields of obj (see note on
%                  skipped field names above). Tables are converted to
%                  structs before hashing.
%               4) none of the above, the hash equals the combined hash of
%                  the hash of each of the properties of obj (see note on
%                  skipped field names above). If obj has no public
%                  properties, the function will fall back to hasing
%                  class(obj)
%  hashStruct - Structure with the hash of each of the fields or properties
%               if obj is either 3 or 4 above. Empty otherwise.
%
%  hash - string of 32 characters with the hexadecimal md5 checksum of
%
% EXAMPLES:
%   obj  = struct('f1', rand(10), 'f2', {{'yes', 'no', struct('a', 2)}});
%   hash = obj2hash(obj);
%
%   G     = computeGeometry(cartGrid([10,10], [100, 100]));
%   rock  = makeRock(G, 100*milli*darcy, 0.1);
%   fluid = initSimpleADIFluid();
%   model = GenericBlackOilModel(G, rock, fluid);
%   [hash, hashStruct] = obj2hash(model);
%
% SEE ALSO:
%   `md5sum`

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
    
    % Optional input arguments
    opt = struct('skip'   , {{}}, ...
                 'maxSize', 1e8 , ...
                 'verbose', []  );
    opt = merge_options(opt, varargin{:});
    % Set verbosity if not sepcified
    if isempty(opt.verbose), opt.verbose = mrstVerbose(); end
    % Check if we are on octave
    octave = mrstPlatform('octave');
    % Function will call itself recursively, so the actual function is
    % implemented in a local version that avoids redoing merge_options,
    % verbosity check, platform check, etc.
    [hash, hashStruct] ...
        = obj2hash_local(obj, opt.skip, opt.maxSize, opt.verbose, octave);
    
end

%-------------------------------------------------------------------------%
function [hash, hashStruct] = obj2hash_local(obj, skip, maxSize, verbose, octave)
% Local version without optional input arguments and checks
    % Initialize second return variable
    hashStruct = [];
    % Avoid empty hash values
    if isempty(obj), obj = '[]'; end
    % Convert tables to structs
    if ~octave && istable(obj), obj = table2struct(obj); end
    % Make function handle for recursive calls
    o2h = @(obj) obj2hash_local(obj, skip, maxSize, verbose, octave);
    if isnumeric(obj) || islogical(obj) || ischar(obj)
        % Object is a type we can hash directly
        if nnz(obj) > maxSize
            if verbose
                warning(['Object has too many nonzero elements (%d), '  , ...
                         'falling back to hashing class name'], nnz(obj));
            end
            obj = class(obj);
        end
        hash = md5sum(obj);
        return
    elseif iscell(obj)
        % Compute hash of each cell array element and hash combined result
        hash = cellfun(@(o) o2h(o), obj, 'UniformOutput', false);
        hash = combineHashes(hash);
        return
    elseif isstruct(obj)
        % Object is a struct. We will have to loop through all fields
        if numel(obj) > 1
            hash = arrayfun(@(o) o2h(o), obj, 'UniformOutput', false);
            hash = combineHashes(hash);
            return
        end
        names = fieldnames(obj);
    else
        % Object is not of a type we can hash directly
        names = properties(obj);
        if isa(obj, 'SimpleTimeStepSelector')
            obj.reset(); % Ensure that timestep selector is reset
        end
    end
    
    if isempty(names)
        % Object has no public properties. Hash the class name instead
        hash = obj2hash(class(obj), 'verbose', false);
        return;
    end
    
    % Exclude any fields we have asked to skip
    names = names(~ismember(names, skip));
    hash  = cell(1, numel(names));
    % Loop through all fields/properties and compute hash
    for i = 1:numel(names)
        % Avoid protected (returned by properties(obj) in Octave)
        if octave && ~isprop(obj, names{i}) && ~isfield(obj, names{i})
            continue
        end
        prop    = obj.(names{i});
        hash{i} = o2h(prop);
    end
    % Extract names that was used in the hash
    keep  = ~cellfun(@isempty, hash);
    names = names(keep);
    hash  = hash(keep);
    % Make struct
    hashStruct = cell2struct(hash, names', 2);
    % Hash the combined result
    hash = combineHashes(hash);
end

%-------------------------------------------------------------------------%
function hash = combineHashes(hashes)
    hash = strjoin(hashes, '_');
    hash = md5sum(hash);
end
