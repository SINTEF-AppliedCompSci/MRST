function tab = preprocessTablePVT(tab)
    nk = numel(tab.key);
    exp = cell(1, nk);
    for tn = 1:nk
        % Precompute each table from the compressed format
        lns = tab.pos(tn):(tab.pos(tn+1)-1);
        T = tab.data(lns, :);
        if size(T, 1) > 1 && T(1, 1) > T(2, 1)
            % If the data is decreasing rather than increasing, we reverse
            % the table.
            T = T(end:-1:1, :);
        end
        exp{tn} = struct('x', T(:, 1), 'F', T(:, 2:end));
    end
    tab.expanded = exp;
    T_sat = tab.data(tab.pos(1:end-1), :);
    tab.sat = struct('x', T_sat(:, 1), 'F', T_sat(:, 2:end));
end