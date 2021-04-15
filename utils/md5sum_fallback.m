function hash = md5sum_fallback(varargin)
    % Alternative implementation of md5sum for systems without C compiler.
    % Requires java.

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

    try
        md = java.security.MessageDigest.getInstance('MD5');
    catch ME
        error('Missing java md5 support!')
    end

    for i = 1:nargin
        addtosum(md, varargin{i});
    end

    digest = md.digest();
    int = java.math.BigInteger(1, digest);
    hash = lower(char(int.toString(16)));
end

function addtosum(md, value)
ui = @(v) typecast(v, 'uint8');

    if issparse(value)
        [i, j, v] = find(value);
        if ~isempty(v)
            md.update(ui(v));
            md.update(ui(i));
            md.update(ui(j));
        end
    elseif isnumeric(value)
        if ~isempty(value)
            md.update(ui(value(:)));
        end
    elseif ischar(value) || islogical(value)
        if ~isempty(value)
            md.update(uint8(value(:)));
        end
    elseif isstruct(value)
        [n m] = size(value);
        for j = 1:m
            for i = 1:n
                f = fieldnames(value(i,j));
                for k = 1:numel(f)
                    addtosum(md, value(i,j).(f{k}));
                end
            end
        end
    elseif iscell(value)
        [n m] = size(value);
        for j = 1:m
            for i = 1:n
                c = value{i,j};
                if ~isempty(c)
                    addtosum(md, value{i});
                end
            end
        end
    else
        warning('Unknown Matlab object.\n')
    end
end
