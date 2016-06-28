function h = plotGridDeformed(G,u,varargin)

    G.nodes.coords = G.nodes.coords + u;
    h = plotGrid(G,varargin{:});

end
