function xv = value(x)
    if isnumeric(x) || islogical(x)
        xv = x;
    elseif iscell(x)
        % Cell arrays are converted to matrices
        x = cellfun(@value, x, 'UniformOutput', false);
        xv = [x{:}];
    elseif isstruct(x)
        fn = fieldnames(x);
        for i = 1:numel(fn)
            f = fn{i};
            x.(f) = value(x.(f));
        end
        xv = x;
    else
        xv = x;
    end
end