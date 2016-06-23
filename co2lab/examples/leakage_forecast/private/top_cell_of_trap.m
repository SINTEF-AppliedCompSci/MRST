function top = top_cell_of_trap(Gt, ta)
% explicitly get top cells of each trap region
% (in case ta.tops is not indexed properly according to trap number)

    num_traps = numel(ta.trap_z);
    top = [];
    for i = 1:num_traps
        tcells = find(ta.traps == i);
        %tz = zeros(Gt.cells.num,1);
        tz = NaN(Gt.cells.num,1);
        tz(tcells) = Gt.cells.z(tcells);
        top(i,1) = find(tz == min(tz)); % tops at at the min elevations of each trap region
        clear tcells
    end


end