function c = getWellPolymerDG(G, W, perf2well, e_w)

    cell2well = nan(G.cells.num,1);
    cell2well(vertcat(W.cells)) = perf2well;
    c = vertcat(W.c);
    c = c(cell2well(e_w));

end