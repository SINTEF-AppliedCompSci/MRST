function graph = getStateFunctionGroupingDependencyGraph(varargin)
% Get dependency graph for one or more StateFunctionGrouping instances
    props = varargin;
    np = numel(props);
    names = cell(np, 1);
    observed = cell(np, 1);
    origin = cell(np, 1);
    group_names = cell(np, 1);
    for i = 1:np
        group = props{i};
        [names{i}, observed{i}] = getNames(group);
        origin{i} = i*ones(numel(names{i}), 1);
        group_names{i} = class(group);
    end
    
    all_names = unique(vertcat(names{:}, observed{:}));
    n = numel(all_names);
    A = zeros(n, n);
    category = zeros(n, 1);
    for i = 1:n
        name =  all_names{i};
        for j = 1:np
            ix = find(strcmp(names{j}, name));
            if ~isempty(ix)
                category(i) = j;
                break
            end
        end
        
        if isempty(ix)
            % Root node
            continue
        end
        A(:, i) = getActive(all_names, props{category(i)}, name);
    end
    graph = struct('C', A,... % Dependency graph
                   'FunctionNames', {all_names}, ... %Names of functions (ordered)
                   'GroupIndex', category, ...% Index of group (with zero for state)
                   'GroupNames', {group_names});% Names of groups (class name)
end

function [names, observed] = getNames(group)
    names = group.getNamesOfStateFunctions();
    observed = {};
    for i = 1:numel(names)
        prop = group.getStateFunction(names{i});
        ext = prop.externals;
        if not(isempty(ext))
            observed = [observed; {ext.name}'];
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

