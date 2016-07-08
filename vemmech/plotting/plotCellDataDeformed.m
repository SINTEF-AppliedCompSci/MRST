function h = plotCellDataDeformed(G, data, u, varargin)
%
%
% SYNOPSIS:
%   function h = plotCellDataDeformed(G, data, u, varargin)
%
% DESCRIPTION: Plot cell data on a grid which is deformed using a given
% displacement field.
%
% PARAMETERS:
%   G        - Grid structure
%   data     - Data to be plotted
%   u        - Displacement field
%   varargin - Optional parameters that are passed further to the function plotCellData
%
% RETURNS:
%   h - handle to plot
%
% EXAMPLE:
%
% SEE ALSO: plotCellData
%

    G.nodes.coords = G.nodes.coords + u;
    if(any(G.faces.areas < 0))
       warning('Deformed grid as negive face areas') 
    end
    if(any(G.cells.volumes < 0))
       warning('Deformed grid has negative volumes')
    end
     
    h = plotCellData(G, data, varargin{:});
end
