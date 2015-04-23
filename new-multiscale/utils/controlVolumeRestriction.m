function restrict = controlVolumeRestriction(partition)
    n = numel(partition);
    m = max(partition);
    i = (1:n) .';
    j = partition(i);
    v = ones(n, 1);
    restrict = sparse(i, j, v, n, m)';
end
