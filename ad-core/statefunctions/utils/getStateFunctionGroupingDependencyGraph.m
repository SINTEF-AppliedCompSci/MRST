function graph = getStateFunctionGroupingDependencyGraph(varargin)
% Get dependency graph for one or more StateFunctionGrouping instances
    props = varargin;
    np = numel(props);
    names = cell(np, 1);
    implementation = cell(np, 1);
    observed = cell(np, 1);
    origin = cell(np, 1);
    group_names = cell(np, 1);
    for i = 1:np
        group = props{i};
        [names{i}, observed{i}, implementation{i}] = getNames(group);
        origin{i} = i*ones(numel(names{i}), 1);
        group_names{i} = class(group);
    end
    all_dependencies = vertcat(observed{:});
    category = rldecode((1:np)', cellfun(@numel, names));
    % External dependencies
    dep = arrayfun(@(x) [x.grouping, '.', x.name], all_dependencies, 'UniformOutput', false);
    [observed_deps, keep_dep] = unique(dep);
    all_dependencies = all_dependencies(keep_dep);
    nobs = numel(observed_deps);
    
    all_names = vertcat(names{:});
    all_implementation = vertcat(implementation{:});
    nnames = numel(all_names);
    n = nnames + nobs;
    A = zeros(n, n);
    for i = 1:nnames
        name = all_names{i};
        c = category(i);
        sf = props{c}.getStateFunction(name);
        % Handle same-group dependencies
        A(1:nnames, i) = ismember(all_names, sf.dependencies);
        % Handle externals
        for j = 1:numel(sf.externals)
            d = sf.externals(j);
            isExt = arrayfun(@(x) strcmp(x.name, d.name) && strcmp(x.grouping, d.grouping), all_dependencies);
            A((nnames+1):end, i) = isExt;
        end
    end
    % Figure out internal and external dependencies
    ext_groups = arrayfun(@(x) x.grouping, all_dependencies, 'UniformOutput', false);
    ext_groups = unique(ext_groups);
    % Number in interior groups - provided here
    n_base_group = numel(group_names);
    % Groups that were induced from dependencies
    extra_categories = setdiff(ext_groups, group_names);
    group_names = [group_names; extra_categories];
    
    nd = numel(all_dependencies);
    dep_names = cell(nd, 1);
    dep_impl = cell(nd, 1);
    dep_cat = zeros(nd, 1);
    for i = 1:nd
        d = all_dependencies(i);
        dep_names{i} = d.name;
        dep_cat(i) = find(strcmp(group_names, d.grouping));
        dep_impl{i} = '?';
    end
    I = [all_implementation; dep_impl];
    N = [all_names; dep_names];
    C = [category; dep_cat];
    % Known dependencies are 1, unknown dependencies are -1
    group_types = 1 - 2*double(C > n_base_group);
    % State is assigned zero
    group_types(strcmp(group_names(C), 'state')) = 0;
    
    assert(numel(I) == numel(C));
    assert(numel(I) == numel(N));
    assert(numel(I) == numel(C));

    combined = cellfun(@(x, y) sprintf('%s.%s', x, y), group_names(C), N, 'UniformOutput', false);
    [~, pos] = uniqueStable(combined);
    
    graph = struct('C', A(pos, pos),... % Dependency graph
                   'Implementation', {I(pos)}, ...% Implementation of specific functions (ordered)
                   'FunctionNames', {N(pos)}, ... % Names of functions (ordered)
                   'GroupIndex', C(pos), ...% Index into GroupNames
                   'GroupTypes', group_types(pos), ... % 1 for connections inside the groups provided, -1 for outside, 0 for state
                   'GroupNames', {group_names});% Names of groups (class name)
end

function [names, observed, impl] = getNames(group)
    if isempty(group)
        error('Empty state function. Did you forget to set it up?');
    end
    [names, ~, impl] = group.getNamesOfStateFunctions();
    observed = [];
    for i = 1:numel(names)
        prop = group.getStateFunction(names{i});
        ext = prop.externals;
        if not(isempty(ext))
            observed = [observed; ext]; %#ok
        end
    end
end

function act = getActive(all_names, group, name)
    property = group.getStateFunction(name);
    dep = ismember(all_names, property.dependencies);
    if not(isempty(property.externals))
        ext = ismember(all_names, {property.externals.name});
    else
        ext = false;
    end
    act = dep | ext;
end

