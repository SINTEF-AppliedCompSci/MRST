function h = plotGridDeformed(G, u, varargin)
%
%
% SYNOPSIS:
%   function h = plotGridDeformed(G, u, varargin)
%
% DESCRIPTION: plot a deformed grid using a given displacement field
%
% PARAMETERS:
%   G        - Grid structure
%   u        - Displacement field
%   varargin - Optional parameters that are passed further to the function plotCellData
%
% RETURNS:
%   h - handle to plot
%
% EXAMPLE:
%
% SEE ALSO:
%

    G.nodes.coords = G.nodes.coords + u;
    h = plotGrid(G, varargin{:});

end
