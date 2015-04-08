function data = coarseDataToFine(CG, data)
% Convert coarse grid dataset into fine grid representation
    if isempty(data)
        return
    end
    if isnumeric(data);
        data = getData(CG, data);
        return
    end
    
    ic = iscell(data);
    
    for i = 1:numel(data)
        if ic
            d = data{i};
        else
            d = data(i);
        end
        
        if isnumeric(d)
            d = getData(CG, d);
        else
            fn = fieldnames(d);
            for j = 1:numel(fn)
                f = fn{j};
                d.(f) = getData(CG, d.(f));
            end
        end
        if ic
            data{i} = d;
        else
            data(i) = d;
        end
    end
end

function data = getData(CG, data)
    sz = size(data);
    if sz(1) == CG.cells.num
        data = data(CG.partition, :);
    elseif all(sz == [1, CG.cells.num])
        data = data(CG.partition);
    end
end
