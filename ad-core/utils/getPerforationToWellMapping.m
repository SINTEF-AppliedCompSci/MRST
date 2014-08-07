function perf2well = getPerforationToWellMapping(w)
    nConn       = cellfun(@numel, {w.cells})'; % # connections of each well
    perf2well   = rldecode((1:numel(w))', nConn);
end