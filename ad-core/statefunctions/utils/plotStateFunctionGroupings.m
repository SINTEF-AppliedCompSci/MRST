function varargout = plotStateFunctionGroupings(props, varargin)
% Plot stateFunctionGrouping
    opt = struct('EdgeColor', [0.3, 0.3, 0.3], ...
                 'LineWidth', 1, ...
                 'Start', {{}}, ...
                 'Stop', {{}}, ...
                 'Center', {{}}, ...
                 'Label', 'name', ...
                 'ColorizeEdges', 'in', ...
                 'Layout', 'layered', ...
                 'includeState', true);
    [opt, arg] = merge_options(opt, varargin{:});
    if ~iscell(props)
        props = {props};
    end
    assert(exist('digraph', 'file') > 0, 'Plotting dependency graphs requires Matlab R2015b or newer.');
    depgraph = getStateFunctionGroupingDependencyGraph(props{:});
    category = depgraph.GroupIndex;
    C = depgraph.C;
    fnames = depgraph.FunctionNames;
    impl = depgraph.Implementation;
    
    names = fnames;
    if opt.includeState
        n = size(C, 1);
        left = zeros(n+1, 1);
        isState = reshape(category == 0, 1, []);
        
        names = ['state'; names];
        impl = ['state'; impl];
        C = [left, [isState; C]];
        category = [-1; category];
        src = 1;
    else
        src = names(category == min(category));
    end
    clear graph
    [G, ~, category, keep] = buildStateFunctionDigraph(C, names, category, ...
                    'Start', opt.Start, 'Stop', opt.Stop, 'Center', opt.Center);
    
    switch lower(opt.Label)
        case 'name'
            % Do nothing
        case 'implementation'
            G.Nodes = impl(keep);
        case 'all'
            local_names = names;
            for i = 1:numel(local_names)
                if ~isempty(impl{i})
                    local_names{i} = sprintf('%s (%s)', local_names{i}, impl{i});
                end
            end
            G.Nodes = local_names(keep);
        otherwise
            error('Bad label %s', opt.Label);
    end
    if strcmpi(opt.Layout, 'layered')
        p = plot(G, 'layout', opt.Layout, 'Sources', src, ...
            'AssignLayers', 'asap', ...
            'direction', 'right', arg{:});
    else
        p = plot(G, 'layout', opt.Layout, arg{:});
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
    
    varargout = cell(1, nargout);
    if nargout > 0
        varargout{1} = p;
        if nargout > 1
            varargout{2} = G;
        end
    end
end
