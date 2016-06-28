function h = plotGridDeformed(G,data,u,varargin)
% plotGridDeformed(G,data,u,varargin)
%{ 
Copyright 2009-2014 SINTEF ICT, Applied Mathematics
%} 
    G.nodes.coords=G.nodes.coords+u;
    %G=computeGeometry(G);
    if(any(G.faces.areas<0))
       warning('Deformed grid as negive face areas') 
    end
    if(any(G.cells.volumes<0))
       warning('Deformed grid has negative volumes')
    end
     
    h= plotCellData(G,data,varargin{:});
end
