function W = combineMSwithRegularWells(W_regular, W_ms)
    % Combine regular and MS wells, accounting for missing fields
    allflds = unique([fieldnames(W_ms); fieldnames(W_regular)]);
    for i = 1:numel(allflds)
        fld = allflds{i};
        for j = 1:numel(W_regular)
            if ~isfield(W_regular(j), fld)
                W_regular(j).(fld) = [];
            end
        end
        for j = 1:numel(W_ms)
            if ~isfield(W_ms(j), fld)
                W_ms(j).(fld) = [];
            end
        end
    end

    W = [W_regular; W_ms];
    for i = 1:numel(W)
        W(i).isMS = i > numel(W_regular);
    end
end