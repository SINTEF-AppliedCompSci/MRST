function hash = struct2hash(opt, varargin)
% Compute hash value of options struct
    % Get field names and corresponding values
    optnames = fieldnames(opt);
    optval   = struct2cell(opt);
    % Keep only fields that are character or numeric arrays
    for i = 1:numel(optval)
        v = optval{i};
        if isempty(v)
            v = '[]';
        elseif ischar(v) || isnumeric(v) || islogical(v)
            v = v(:)';
        elseif isstruct(v)
            v = arrayfun(@(v) struct2hash(v), v, 'UniformOutput', false);
            v = strjoin(v, '_');
        else
            v = obj2hash(v);
        end
        optval{i} = v;
    end
    % Convert numeric arrays to strings
    fix         = cellfun(@(v) ~ischar(v), optval);
    optval(fix) = cellfun(@(v) num2str(v), optval(fix), ...
                                    'UniformOutput', false);
    % Indices and values
    ix  = [1:2:2*numel(optnames), 2:2:2*numel(optnames)];
    rhs = [optnames; optval];
    % Prepend with varargin if given
    if nargin > 1
        assert(numel(varargin) == 1);
        ix  = [1, ix + 1];
        rhs = [varargin{1}; rhs];
    end
    % Make string
    str     = cell(1, 2*numel(optval) + (nargin>1));
    str(ix) = rhs;
    str     = horzcat(str{:});
    % Compute hash value of string
    hash = str2hash(str);
end