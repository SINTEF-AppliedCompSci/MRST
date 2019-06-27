function varargout = plotStateFunctionGroupings(props, varargin)
% Plot stateFunctionGrouping
    opt = struct('EdgeColor', [0.3, 0.3, 0.3], ...
                 'LineWidth', 1, ...
                 'Start', {{}}, ...
                 'Stop', {{}}, ...
                 'Center', {{}}, ...
                 'Layout', 'layered', ...
                 'includeState', true);
    [opt, arg] = merge_options(opt, varargin{:});
    if ~iscell(props)
        props = {props};
    end
    assert(exist('digraph', 'file') > 0, 'Plotting dependency graphs requires Matlab R2015b or newer.');
    graph = getStateFunctionGroupingDependencyGraph(props{:});
    category = graph.GroupIndex;
    C = graph.C;
    names = graph.FunctionNames;
    if opt.includeState
        n = size(C, 1);
        left = zeros(n+1, 1);
        isState = reshape(category == 0, 1, []);
        
        names = ['state'; names];
        C = [left, [isState; C]];
        category = [-1; category];
        src = 1;
    else
        src = names(category == min(category));
    end
    clear graph
    [G, names, category] = buildStateFunctionDigraph(C, names, category, ...
                    'Start', opt.Start, 'Stop', opt.Stop, 'Center', opt.Center);
    
    
    if strcmpi(opt.Layout, 'layered')
        p = plot(G, 'layout', opt.Layout, 'Sources', src, arg{:});
    else
        p = plot(G, 'layout', opt.Layout, arg{:});
    end
    p.NodeCData = category;
    p.LineWidth = opt.LineWidth;
    p.EdgeColor = opt.EdgeColor;
    colormap(lines(max(category) - min(category) + 1))
    
    if nargout
        varargout = {p};
    end
end
