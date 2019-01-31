function xv = value(x)
    if isnumeric(x) || islogical(x)
        xv = x;
    elseif iscell(x)
        sz = size(x);
        xv = cellfun(@value, x, 'UniformOutput', false);
        if sz(1) == 1 && sz(2) > 1
            % Cell arrays are converted to matrices
            xv = [xv{:}];
        end
    elseif isstruct(x)
        fn = fieldnames(x);
        for i = 1:numel(fn)
            f = fn{i};
            for j = 1:numel(x)
                x(j).(f) = value(x(j).(f));
            end
        end
        xv = x;
    else
        xv = x;
    end
end