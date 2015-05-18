function printConvergenceReport(names, values, converged, iteration)
    nl = cellfun(@numel, names);
    sep = repmat('=', 1, sum(nl) + 3*numel(nl) + 1);
    if iteration == 1
        fprintf('%s\n', sep);
        fprintf('| %-8s ', names{:});
        fprintf('|\n%s\n', sep);
    end
    for i = 1:numel(values)
        linen = max(nl(i), 8);
        fprintf(['| %-', num2str(linen), '.2e '], values(i));
    end
    fprintf('|\n')
    if converged
        fprintf('%s\n', sep);
    end
end
