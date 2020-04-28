function bnds = processBounds(W, varargin)
nw   = numel(W);
flds = {'bhp', 'wrat', 'orat', 'grat', 'lrat', 'rate'};
inp  = [flds; repmat({nan(nw,1)}, [1, numel(flds)]) ];
[lower, upper] = deal(struct(inp{:}));
[inj, injectors]  = deal(vertcat(W.sign)>0); %#ok
[prod, producers] = deal(vertcat(W.sign)<0); %#ok
for k =1:2:numel(varargin)
    [expr , vals] = deal(varargin{k:k+1});
    assert(isa(expr, 'char') && isnumeric(vals));
    try
        eval(['lower.', expr, ';']);
    catch
        error('Unknown expression: %s ', expr);
    end
    eval(['lower.', expr, ' = vals(:,1);']);
    eval(['upper.', expr, ' = vals(:,2);']);
end
bnds = lower;
for k = 1:numel(flds)
    bnds.(flds{k}) = [lower.(flds{k}), upper.(flds{k})];
end
% if producer bounds are positive, assume absolute, so change sign and pos 
pflds =  {'wrat', 'orat', 'grat', 'lrat'};
for k = 1:numel(pflds)
    ix = any(bnds.(pflds{k})>0, 2);
    if any(ix)
        bnds.(pflds{k})(ix,:) = -bnds.(pflds{k})(ix, [2 1]);
    end
end     
end