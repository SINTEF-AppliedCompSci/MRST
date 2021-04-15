function graph = getStateFunctionGroupingDependencyGraph(varargin)
% Get dependency graph for one or more StateFunctionGrouping instances

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    function_groupings = varargin;
    np = numel(function_groupings);
    names = cell(np, 1);
    implementation = cell(np, 1);
    observed = cell(np, 1);
    origin = cell(np, 1);
    group_names = cell(np, 1);
    symbols = cell(np, 1);
    for i = 1:np
        group = function_groupings{i};
        [names{i}, observed{i}, implementation{i}, symbols{i}] = getNames(group);
        origin{i} = i*ones(numel(names{i}), 1);
        group_names{i} = class(group);
    end
    category = rldecode((1:np)', cellfun(@numel, names));
    all_dependencies = vertcat(observed{:});
    all_names = vertcat(names{:});
    all_symbols = vertcat(symbols{:});
    all_implementation = vertcat(implementation{:});
    names_full = cellfun(@(x, y) sprintf('%s.%s', x, y),...
                    group_names(category), all_names, 'UniformOutput', false);
    % Strip away subclass relationships. The dependencies are documented
    % via base classes while we may be dealing with a superclass.
    all_dependencies = fixSubclassDependencies(function_groupings, all_dependencies);
    clear observed
    % External dependencies
    dep_full = arrayfun(@(x) [x.grouping, '.', x.name], all_dependencies, 'UniformOutput', false);
    [observed_deps, keep_dep] = setdiff(dep_full, names_full);

    all_dependencies = all_dependencies(keep_dep);
    nobs = numel(observed_deps);
    
    nnames = numel(all_names);
    n = nnames + nobs;
    A = zeros(n, n);
    for i = 1:nnames
        name = all_names{i};
        c = category(i);
        sf = function_groupings{c}.getStateFunction(name);
        
        % Handle same-group dependencies
        local = category == c;
        A(1:nnames, i) = ismember(all_names, sf.dependencies) & local;
        % Handle externals
        externals = fixSubclassDependencies(function_groupings, sf.externals);
        for j = 1:numel(sf.externals)
            d = externals(j);
            % Match to dependencies we have not observed
            isExt = arrayfun(@(x) strcmp(x.name, d.name) && strcmp(x.grouping, d.grouping), all_dependencies);
            A((nnames+1):end, i) = A((nnames+1):end, i) | isExt;
            % Match to dependencies which are known as a part of the set of
            % state function groupings provided for this function
            groupMatch = strcmp(group_names, d.grouping);
            if any(groupMatch)
                isNor = ismember(all_names, d.name) & category == find(groupMatch);
                A(1:nnames, i) = A(1:nnames, i) | isNor;
            end
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
    dep_symbl = cell(nd, 1);
    dep_cat = zeros(nd, 1);
    for i = 1:nd
        d = all_dependencies(i);
        dep_names{i} = d.name;
        dep_cat(i) = find(strcmp(group_names, d.grouping));
        if strcmpi(d.grouping, 'state')
            dep_impl{i} = 'state';
        else
            dep_impl{i} = '?';
        end
        dep_symbl{i} = d.name;
    end
    I = [all_implementation; dep_impl];
    N = [all_names; dep_names];
    S = [all_symbols; dep_symbl];
    C = [category; dep_cat];
    % Known dependencies are 1, unknown dependencies are -1
    group_types = 1 - 2*double(C > n_base_group);
    % State is assigned zero
    group_types(strcmp(group_names(C), 'state')) = 0;
    
    assert(numel(I) == numel(C));
    assert(numel(I) == numel(N));
    assert(numel(I) == numel(S));

    graph = struct('C', A,                     ... % Dependency graph
                   'Implementation', {I},      ... % Implementation of specific functions (ordered)
                   'FunctionNames', {N},       ... % Names of functions (ordered)
                   'Labels',       {S},        ... % Label/shortname for each element
                   'GroupIndex', C,            ... % Index into GroupNames
                   'GroupTypes', group_types,  ... % 1 for connections inside the groups provided, -1 for outside, 0 for state
                   'GroupNames', {group_names});   % Names of groups (class name)
end

function deps = fixSubclassDependencies(groups, deps)
    for i = 1:numel(deps)
        for j = 1:numel(groups)
            if isa(groups{j}, deps(i).grouping)
                deps(i).grouping = class(groups{j});
                break
            end
        end
    end
end

function [names, observed, impl, symbols] = getNames(group)
    if isempty(group)
        error('Empty state function. Did you forget to set it up?');
    end
    [names, ~, impl, symbols] = group.getNamesOfStateFunctions();
    observed = [];
    for i = 1:numel(names)
        prop = group.getStateFunction(names{i});
        ext = prop.externals;
        if not(isempty(ext))
            observed = [observed; ext]; %#ok
        end
    end
end
