function [groupSols, wellSols] = computeWellGroupSols(wellSols, W, varargin)
    
    groups = getGroups(W);

    n = numel(wellSols);
    groupSols = cell(n,1);
    for i = 1:n
        wellSol  = wellSols{i};
        groupSol = aggregateGroupData(groups, wellSol);
        groupSols{i} = groupSol;
        wellSols{i}  = [wellSol, groupSol];
    end
    
end

%-------------------------------------------------------------------------%
function groups = getGroups(W)
    
    names = {W.group};
    groups = struct('name', unique(names));
    for i = 1:numel(groups)
        groups(i).mask = strcmpi(names, groups(i).name);
    end

end

%-------------------------------------------------------------------------%
function groupSol = aggregateGroupData(groups, wellSol)

    sample = wellSol(1);
    groupSol = groups;
    
    rateNames = {'qWs', 'qOs', 'qGs'};
    active    = ~cellfun(@isempty, {sample.qWs, sample.qOs, sample.qGs});
    dummy      = rateNames(~active);
    rateNames = rateNames(active);
    rateNames{end+1} = 'effect';
    
    dummy = [dummy, {'type', 'val'}];
    skip = {'name'};
    
    fnames = reshape(fieldnames(sample), 1, []);
    for i = 1:numel(groups)
        wsg = wellSol(groups(i).mask);
        for fname = fnames
            if any(strcmpi(fname{1}, skip))
                continue
            elseif any(strcmpi(fname{1}, dummy)) || ...
                    numel(wsg(1).(fname{1})) ~= 1
                    groupSol(i).(fname{1}) = wsg(1).(fname{1});
                continue
            elseif any(strcmpi(fname{1}, rateNames))
                get = @(v) sum(v);
            elseif strcmpi(fname{1}, 'status')
                get = @(v) nnz(v)./numel(v);
            elseif strcmpi(fname{1}, 'bhp')
                get = @(v) mean(v);
            else
                get = @(v) v(1);
            end
            v = [wsg.(fname{1})];
            groupSol(i).(fname{1}) = get(v);
        end
    end
    
    if isfield(wellSol(1), 'T')
        for i = 1:numel(groups)
            wsg = wellSol(groups(i).mask);
            q   = [wsg.qWs];
            T   = [wsg.T];
            T   = sum(q.*T)./sum(q);
            groupSol(i).T = T;
        end
    end

    groupSol = arrayfun(@(gs) rmfield(gs, 'mask'), groupSol);
    
end