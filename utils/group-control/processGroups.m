function groups = processGroups(model, groups, state, pot)

    ng         = numel(groups);
    groups0    = groups;
    groupNames = cellfun(@(group) group.name, groups0, 'UniformOutput', false);
    isChild = false(1, ng);
    for i = 1:ng
        if isChild(i), continue; end
        group = groups0{i};
        groups{i} = {group};
        if ~isfield(group, 'children'), continue; end
        if strcmpi(group.type, 'none'), continue; end
        children = group.children;
        nc = numel(children);
        childrenPot = zeros(1, nc);
        active = true(nc,1);
        for j = 1:nc
            child = children{j};
            mask = getGroupMask(model, state, child);
            if ~any(mask), active(j) = false; continue; end
            chix = strcmpi(groupNames, child);
            if any(chix)
                children{j} = groups0{chix};
                children{j}.type = group.type;
                children{j}.val  = group.val;
                children{j}.T    = group.T;
                isChild(chix)    = true;
            else
                children{j} = struct('name', child     , ...
                                     'type', group.type, ...
                                     'val' , group.val , ...
                                     'frac', nan       , ...
                                     'T'   , group.T   );
            end
            childrenPot(j) = sum(pot(mask));
        end
        groupPot = sum(childrenPot);
        for j = 1:nc
            if ~active(j), continue; end
            if ~isfield(children{j}, 'frac') || isnan(children{j}.frac)
                children{j}.frac = childrenPot(j)./groupPot;
            end
        end
        children = children(active);
        if any(active)
            groups{i} = children;
        else
            groups{i}{1}.type = 'none';
        end
    end
    
    groups = groups(~isChild);
    groups = horzcat(groups{:});

    ng = numel(groups);
    for i = 1:ng
        group = groups{i};
        if strcmpi(group.type, 'none'), continue; end
        if ischar(group.val)
            gix = strcmpi(group.val, groupNames);
            if isfield(groups0{gix}, 'children')
                group.val = groups0{gix}.children;
            end
        end
        groups{i} = group;
    end
    
end