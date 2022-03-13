function [s, h, hmax]  = normalizeValuesVE(g_top, sol, fluid, varargin)
% Normalizes values for different VE grids and formulations.
%
%
% SYNOPSIS:
%   [s h hmax]  = normalizeValuesVE(g_top, sol, fluid)
%   [s h hmax]  = normalizeValuesVE(g_top, sol, fluid, 'CoupledGrid', g_coupled)
%
% PARAMETERS:
%   g_top   - top surface 2D grid as defined by function 'topSurfaceGrid'.
%   
%   sol     - Reservoir solution from for example initResSol
%
%   fluid   - A valid fluid object for the formulation
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.
%             Supported options are:
%               'CoupledGrid' -- The coupled grid. Needed if using a
%               coupled formulation.
%               'Rock' -- fine scale rock structure for coupled h
%               calculations using sat2height
%
% RETURNS:
%   s       - Fine scale saturation based on the values. This will be
%             suitable for plotting with the full 3D grid the top surface
%             grid was made from. 
%
%   h       - Height values suitable for plotting with the top surface
%             grid. Contains g_top.cells.num elements.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

   opt = struct('CoupledGrid', [], 'Rock', []);
   opt = merge_options(opt, varargin{:});
   
   if ~isempty(opt.CoupledGrid)
        g_coupled = opt.CoupledGrid;
        % Find the h values everywhere in coupled grid
        h_tot = fluid.sat2height(sol);
        
        h = zeros(g_top.cells.num,1);
        % Find indices in the top surface grid corresponding to the columns
        % of g_coupled where VE is used
        colind = g_coupled.cells.mapTopSurface>0; %& ~g_coupled.cells.region3D;
        
        % Indices of 3d cells in the coupled grid in the full 3d grid
        full3dind = g_coupled.cells.inx3D(g_coupled.cells.region3D);
        
        % Map the top values to columns so that we return the correct
        % values
        h(colind) = h_tot(colind);
        
        tmp.h = h;
        
        tmp.h_max = h; %TODO fiks henting av historikk fra faktisk sol
        % Find saturation based on the constructed height values
        s = height2Sat(tmp, g_top, fluid); 
        % However, in the parts of the grid where we already have a fine
        % scale saturation we should use what we already got.
        s(full3dind) = sol.s(g_coupled.cells.region3D);
        if ~isempty(opt.Rock)
           h = sat2height(s, g_top, opt.Rock);
        end
        
   elseif isfield(fluid, 'sat2height')
        % We probably have a S fluid
        % Use top surface saturation to find height 
        [h, hmax] = fluid.sat2height(sol);
        sol.h = h;
        % Fine scale saturation
        s = height2Sat(sol, g_top, fluid);
    elseif any(strcmpi(g_top.type, 'topsurfacegrid'))
        % We have a H fluid
        h = sol.h;
        hmax = sol.h_max;
        s = height2Sat(sol, g_top, fluid);
    end
    assert(numel(s)==numel(g_top.columns.cells));
    assert(numel(h)==g_top.cells.num);
end
