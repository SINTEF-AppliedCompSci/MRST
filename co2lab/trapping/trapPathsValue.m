function values = trapPathsValue(G, trees, res)
    alreadyPlotted = false(G.cells.num, 1);
    values = zeros(G.cells.num, 1);
    for i = 1:numel(trees)
        cl = res.cell_lines(trees(i).traps);
        lines = false(G.cells.num, 1);
        for j = 1:numel(cl)
            tmp = cl{j};
            if ~isempty(tmp)
                lines(tmp{1}) = true;
            end
        end
        
        c = (lines | ismember(res.traps, trees(i).traps)) ... 
            & ~alreadyPlotted;
        
        alreadyPlotted(c) = true;
        values(c) = trees(i).value;
    end
end
