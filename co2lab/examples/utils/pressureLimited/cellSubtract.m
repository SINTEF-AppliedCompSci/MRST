function c = cellSubtract(a, b)
    % we assume that a and b should be single-indexed
    c = a;
    for i = 1:numel(c)
        c{i} = c{i} - b{i};
    end
end