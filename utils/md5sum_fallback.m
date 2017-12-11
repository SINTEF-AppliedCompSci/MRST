function hash = md5sum_fallback(varargin)
    % Alternative implementation of md5sum for systems without C compiler.
    % Requires java.
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
