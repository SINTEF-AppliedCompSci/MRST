function h = plotCellDataDeformed(G, data, u, varargin)

    G.nodes.coords = G.nodes.coords + u;
    if(any(G.faces.areas < 0))
       warning('Deformed grid as negive face areas') 
    end
    if(any(G.cells.volumes < 0))
       warning('Deformed grid has negative volumes')
    end
     
    h = plotCellData(G, data, varargin{:});
end
