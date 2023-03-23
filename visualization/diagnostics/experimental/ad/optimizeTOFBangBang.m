function [D_best W_best history objlist] = optimizeTOFBangBang(optimizeHandle, minvals, W, targets)
    [D_best W_best history] = optimizeHandle(W);

    if numel(minvals) == 1
        minvals = minvals*ones(numel(targets), 1);
    end

    totSum = sum([W(targets).val]);

    objlist = history.gradient(end).objective.val;
    for i = 1:numel(targets)
        noti = setdiff(1:numel(targets), i);
        W(i).val = totSum - sum(minvals(targets(noti)));
        for j = noti
            W(j).val = minvals(targets(j));
        end
        [D W h] = optimizeHandle(W);
        if history.gradient(end).objective.val > h.gradient(end).objective.val
            D_best = D;
            W_best = W;
            history = h;
        end
        objlist = [objlist; h.gradient(end).objective.val];
    end
end
