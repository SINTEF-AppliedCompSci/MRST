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
% SEE ALSO: `plotCellData`
%

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

    G.nodes.coords = G.nodes.coords + u;
    if(any(G.faces.areas < 0))
       warning('Deformed grid as negive face areas') 
    end
    if(any(G.cells.volumes < 0))
       warning('Deformed grid has negative volumes')
    end
     
    h = plotCellData(G, data, varargin{:});
end
