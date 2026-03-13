function W = addDefaultWellFields(W)
for k =1:numel(W)
    np = numel(W(k).cells);
    W(k).topo = [(0:(np-1))', (1:np)'];
    W(k).status = true;
    W(k).cstatus = true(np, 1);
    W(k).lims = [];
end
