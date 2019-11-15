function v = expandIfUniform(v)
    % Utility which reverses "value" compaction. If given a matrix (logical
    % or numerical) as input, it will expand it to a cell array of vectors
    % such that value(expandIfUniform(x)) is equal to x.
    if (isnumeric(v) || islogical(v)) && size(v, 2) > 1
        n = size(v, 2);
        out = cell(1, n);
        for i = 1:n
            out{i} = v(:, i);
        end
        v = out;
    end
end
