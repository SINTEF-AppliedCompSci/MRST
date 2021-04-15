function varargout = plotStateFunctionGroupings(props, varargin)
% Plot stateFunctionGrouping

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

    opt = struct('EdgeColor', [0.3, 0.3, 0.3], ...
                 'LineWidth', 1, ...
                 'Start', {{}}, ...
                 'Stop', {{}}, ...
                 'Center', {{}}, ...
                 'TextArg', {{}}, ...
                 'Label', 'name', ...
                 'Sources', [], ...
                 'ColorizeEdges', 'in', ...
                 'Layout', 'layered', ...
                 'includeState', true);
    [opt, arg] = merge_options(opt, varargin{:});
    if isa(props, 'PhysicalModel')
        % We got a model! Validate it and output the groupings
        props = props.validateModel();
        props = props.getStateFunctionGroupings();
    end
    if ~iscell(props)
        props = {props};
    end
    assert(exist('digraph', 'file') > 0, 'Plotting dependency graphs requires Matlab R2015b or newer.');
    depgraph = getStateFunctionGroupingDependencyGraph(props{:});
    category = depgraph.GroupIndex;
    C = depgraph.C;
    names = depgraph.FunctionNames;
    groupnames = depgraph.GroupNames;
    impl = depgraph.Implementation;
    labl = depgraph.Labels;
    full_names = cellfun(@(x, y) sprintf('%s.%s', x, y), groupnames(category), names, 'UniformOutput', false); 
    if  ~isempty(opt.Sources)
        src = opt.Sources;
        for i = 1:numel(src)
            if ~contains(src{i}, '.')
                tmp = find(contains(full_names, src{i}), 1, 'first');
                src{i} = full_names{tmp}; % Hope that we match
            end
        end
    elseif opt.includeState
        n = size(C, 1);
        left = zeros(n+1, 1);
        isState = reshape(depgraph.GroupTypes <= 0, 1, []);
        
        names = ['state'; names];
        impl = ['state'; impl];
        labl = ['state'; labl];
        full_names = ['state'; full_names];
        C = [left, [isState; C]];
        category = [1; category+1];
        src = 1;
    else
        types = depgraph.GroupTypes;
        if max(types) == min(types)
            src = full_names{1};
        else
            src = full_names(types == min(types));
        end
    end
    clear graph
    [G, ~, category, keep] = buildStateFunctionDigraph(C, full_names, category, ...
                    'Start', opt.Start, 'Stop', opt.Stop, 'Center', opt.Center, 'FilterNames', names);
    if iscell(src)
        src = setdiff(src, full_names(~keep));
        if isempty(src)
            src = 1;
        end
    end
    % Set layout
    plot_defaults = {'EdgeAlpha', 1};
    switch lower(opt.Layout)
        case 'layered'
            p = plot(G, 'layout', opt.Layout, 'Sources', src, ...
                'AssignLayers', 'asap', ...
                'direction', 'right', plot_defaults{:}, arg{:});
        case 'mrst'
            p = plot(G, 'layout', 'layered', plot_defaults{:}, arg{:});
            p.XData = category;
            d = ones(numel(category), 1);
            cts = accumarray(category, 1);
            for i = 1:max(category)
                local = category == i;
                d(local) = linspace(1, max(cts), sum(local))';
                dx = (1:sum(local))/sum(local);
                p.XData(local) = p.XData(local) + 2*dx;
            end
            p.YData = d;
        otherwise
            p = plot(G, 'layout', opt.Layout, plot_defaults{:}, arg{:});
    end
    % Set labels
    switch lower(opt.Label)
        case 'name'
            p.NodeLabel = names(keep);
        case 'label'
            p.NodeLabel = labl(keep);
        case 'implementation'
            p.NodeLabel = impl(keep);
        case 'all'
            local_names = names;
            for i = 1:numel(local_names)
                if strcmpi(local_names{i}, 'state')
                    % Do nothing
                elseif ~isempty(impl{i})
                    local_names{i} = sprintf('%s (%s)', local_names{i}, impl{i});
                end
            end
            p.NodeLabel = local_names(keep);
        case 'debug'
            % Do nothing, show actual node names
        otherwise
            error('Bad label %s', opt.Label);
    end

    ce = lower(opt.ColorizeEdges);
    switch ce
        case {'in', 'out', 'avg'}
            w = category;
            switch ce
                case 'in'
                    w_fn = @(in, out) w(in);
                case 'out'
                    w_fn = @(in, out) w(out);
                case 'avg'
                    w_fn = @(in, out) 0.5*(w(in) + w(out));
            end
            edges = G.Edges;
            ne = size(edges, 1);
            cData = zeros(ne, 1);
            for edgeNo = 1:ne
                e = edges{edgeNo, 1};
                in = strcmp(G.Nodes.Name, e{1});
                out = strcmp(G.Nodes.Name, e{2});
                cData(edgeNo) = w_fn(in, out);
            end
            set(p, 'EdgeCData', cData);
        case 'none'
            p.EdgeColor = opt.EdgeColor;
        otherwise
            error('Unknown ColorizeEdges option %s. Valid choices: in, out, avg, none', ce);
    end
    
    p.NodeCData = category;
    p.LineWidth = opt.LineWidth;
    
    colormap(lines(max(category) - min(category) + 1))
    if ~isempty(opt.TextArg)
        % Replace labels with custom text
        x = get(p, 'XData');
        y = get(p, 'YData');
        labels = get(p, 'NodeLabel');
        th = text(x, y, labels);
        for j = 1:numel(th)
            set(th(j), opt.TextArg{:})
        end
        % Remove existing labels
        set(p, 'NodeLabel', {});
        % Store text handles in UserData
        set(p, 'UserData', th);
    end
    varargout = cell(1, nargout);
    if nargout > 0
        varargout{1} = p;
        if nargout > 1
            varargout{2} = G;
            if nargout > 2
                varargout{3} = depgraph;
            end
        end
    end
end
