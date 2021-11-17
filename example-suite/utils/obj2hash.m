function hash = obj2hash(obj, varargin)
    if isstruct(obj)
        str = cell(1, numel(obj));
        for i = 1:numel(obj)
            str{i} = struct2hash(obj(i));
        end
        str = strjoin(str, '_');
    elseif ischar(obj)
        str = obj;
    elseif isempty(obj)
        str = '[]';
    elseif isnumeric(obj) || islogical(obj)
        str = num2str(obj(:)');
    elseif isa(obj, 'function_handle')
        str = func2str(obj);
    else
        props = properties(obj)';
        np = numel(props);
        str    = cell(1, np+1);
        str{1} = class(obj);
        for i = 1:np
            str{i+1} = obj2hash(obj.(props{i}));
        end
        str = strjoin(str, '_');
    end
    hash = str2hash(str);
end