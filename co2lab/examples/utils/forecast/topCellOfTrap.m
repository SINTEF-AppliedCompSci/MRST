function tops = topCellOfTrap(Gt, ta)
% explicitly get top cells of each trap region
% (in case ta.tops is not indexed properly according to trap number)

    num_traps = numel(ta.trap_z);
    tops = [];
    for i = 1:num_traps
        tcells = find(ta.traps == i);
        %tz = zeros(Gt.cells.num,1);
        tz = NaN(Gt.cells.num,1);
        tz(tcells) = Gt.cells.z(tcells);
        top = find(tz == min(tz)); % top at the min elevation of trap i
        tops(i,1) = top(1);  % the top elevation may occur at more than one cell, so we only
                             % need to take one cell to get the top
        clear tcells
    end


end