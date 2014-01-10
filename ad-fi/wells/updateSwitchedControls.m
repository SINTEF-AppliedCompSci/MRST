function W = updateSwitchedControls(sol, W)
% Check if producers are becoming injectors and vice versa
% Nothing is done here for now, these should be shut down
wsg = vertcat(W(:).sign);
ssg = vertcat(sol(:).sign);
inx = find(wsg~=ssg);
for k = 1:numel(inx)
    tps  = {'injector', 'producer'};
    tpIx = [1 2];
    if wsg(inx(k)) < 1, tpIx = [2 1];end

    warning(['Well ', W(inx(k)).name, ' has switched from ',tps(tpIx(1)), ' to ' tps(tpIx(2)),'.']);
end

% Check if well-controls have been switch, if so, update W
inx = find(~arrayfun(@(x,y)strcmp(x.type,y.type), W(:), sol(:)));
for k = 1:numel(inx)
    fromTp = W(inx(k)).type;
    toTp   = sol(inx(k)).type;
    fprintf(['Well ', W(inx(k)).name, ' has switched from ', fromTp, ' to ', toTp, '.\n']);
    W(inx(k)).type = toTp;
    W(inx(k)).val  = sol(inx(k)).val;
end
end




