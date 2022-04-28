function groups = processGroups(model, groups, state, pot)

    ng = numel(groups);
    groups0 = groups;
    groupNames = cellfun(@(group) group.name, groups0, 'UniformOutput', false);
    isChild = false(1, ng);
    keep = true(ng,1);
    for i = 1:ng
        if isChild(i), continue; end
        group = groups0{i};
        groups{i} = {group};
        if ~isfield(group, 'children'), continue; end
        if strcmpi(group.type, 'none'), continue; end
        children = group.children;
        nsg = numel(children);
        childrenPot = zeros(1, nsg);
        
        for j = 1:nsg
            child = children{j};
            chix = strcmpi(groupNames, child);
            if any(chix)
                children{j} = groups{chix};
                children{j}.val = group.val;
                children{j}.T   = group.T;
                isChild(chix) = true;
            else
                mask = getGroupMask(model, state, child);
                childrenPot(j) = sum(pot(mask));
                children{j} = struct('name', child      , ...
                    'type'   , group.type, ...
                    'val'    , group.val , ...
                    'frac', nan           , ...
                    'T'      , group.T       );
            end
        end
        groupPot = sum(childrenPot);
        for j = 1:nsg
            if isnan(children{j}.frac)
                children{j}.frac = childrenPot(j)./groupPot;
            end
        end
        groups{i} = children;
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