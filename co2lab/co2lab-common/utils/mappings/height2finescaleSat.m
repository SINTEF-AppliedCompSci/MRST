function [s, smax, seff] = height2finescaleSat(h, hmax, Gt, sw, sg, varargin)
% Convert from height to fine scale saturation.  By default, a sharp-interface
% approximation is assumed.  If a capillary fringe model should be used, the
% function needs to be provided with the optional arguments 'invPc3D', 'rhoW'
% and 'rhoG'.
% 
% SYNOPSIS:
%   s = height2finescaleSat(h, hmax, Gt, sw, wg, varargin)
%
% PARAMETERS:
%   h - CO2 plume thickness.  One scalar value for each column in the
%       top-surface grid.
%
%       Values less than zero are treated as zero while values below the
%       bottom of a column are treated as the column depth.
%
%   hmax - historically maximum thickness.  One scalar value for each
%          column in the top surface grid
%    
%   Gt - A top-surface grid as defined by function 'topSurfaceGrid'.
%
%   sw - residual water saturation
%   sg - residual gas saturation
% 
%   invPc3D (optional) - If this argument is provided, a capillary fringe
%                        model is assumed.  'invPc3D' should then be the
%                        inverse fine-scale capillary pressure function,
%                        i.e. taking a capillary pressure value as argument
%                        and returning the corresponding saturation.  This
%                        function will be available from the fluid model if
%                        a capillary fringe model is used.
%
%   rhoW, rhoG (optional) - water and CO2 densities, given for each cell in Gt.
%                           These are only required input if a capillary fringe
%                           model is called for, i.e. 'invPc3D' provided.
%
% RETURNS:
%   s - Fine scale saturation, one value for each cell in the underlying 3D model.
%       Corresponds to state.s for the 3D problem.
%   smax - Maximum historical value of fine-scale saturation, one value per
%          3D cell.  smax >= s.
%   seff - "Effective" fine scale saturation, i.e. the part of fine-scale
%          saturation that is not considered 'residually trapped', and
%          contributes to flow.  seff <= s <= smax
%
% SEE ALSO:
%   `accumulateVertically`, `integrateVertically`

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = merge_options(struct('invPc3D', [], 'rhoW', [], 'rhoG', []), varargin{:});
    
    if isempty(opt.invPc3D)
        [s, smax, seff] = sharp_interface_h2s(h, hmax, Gt, sw, sg);
    else
        assert(~isempty(opt.rhoW) && ~isempty(opt.rhoG));
        [s, smax, seff] = cap_fringe_h2s(h, hmax, Gt, sw, sg, opt.invPc3D, opt.rhoW, opt.rhoG);
    end
end

% ----------------------------------------------------------------------------
function result = remap(vals, ixs)
    result = zeros(numel(vals), 1);
    result(ixs) = vals;
end

% ----------------------------------------------------------------------------
function [s, smax, seff] = cap_fringe_h2s(h, hmax, Gt, sw, sg, invPc3D, rhoW, rhoG)
    
    % endpoint scaling factor
    C = sg / (1 - sw);
    
    % remap upscaled variables to fine-scale grid
    %remap = @(x, ixs) x(ixs);
    to_finescale = @(var) remap(rldecode(var, diff(Gt.cells.columnPos)), Gt.columns.cells);
    
    iface_depth_all = to_finescale(h + Gt.cells.z); % depth of interface in the column
    iface_depth_max_all = to_finescale(hmax + Gt.cells.z);
    drho_all = to_finescale(rhoW - rhoG);
    
    celltops    = Gt.parent.cells.centroids(:,3) - remap(Gt.columns.dz, Gt.columns.cells) / 2; 
    cellbottoms = Gt.parent.cells.centroids(:,3) + remap(Gt.columns.dz, Gt.columns.cells) / 2; 
    
    % compute capillary pressure and take inverse to get effective and max saturations
    seff = compute_cell_saturations(celltops, cellbottoms, iface_depth_all, drho_all, invPc3D);
    smax = compute_cell_saturations(celltops, cellbottoms, iface_depth_max_all, drho_all, invPc3D);
    
    % combine seff and smax to get current fine-scale saturation
    s = (1 - C) * seff + C * smax;
end

% ----------------------------------------------------------------------------
function s = compute_cell_saturations(celltops, cellbottoms, iface_depth, drho, invPc3D)
    % Compute cell saturation as an average of the cell saturations at a number of
    % internal vertical levels inside the cell.
        
    N = 10; % internal subdivisions for reconstruction

    s = zeros(numel(celltops), 1);
    
    for i = 1:N
        w = 1;
        if i == 1 || i == N
            w = 0.5; % simpson's rule
        end
        level = (celltops * (N-i) + cellbottoms * (i-1)) / (N-1); 
        s_level = 1 - invPc3D(max(iface_depth - level, 0) .* drho .* norm(gravity));
        
        s = s + w * s_level;
    end
    s = s / (N-1);
end

% ----------------------------------------------------------------------------
function [s, smax, seff] = sharp_interface_h2s(h, hmax, Gt, sw, sg)
    
    [s, seff] = deal(zeros(numel(Gt.columns.cells),1));
    
    % n: number of completely filled cells
    % t: fill degree for the column's single partially filled cell
    [n, t] = fillDegree(h, Gt); %
    
    % number of cells in each column
    nc = diff(Gt.cells.columnPos);
    
    % compute internal cellNumber in the column for each cell
    cellNoInCol = mcolon(ones(Gt.cells.num,1), diff(Gt.cells.columnPos))';
    
    % f(cells with s == 1)    > 0
    % f(cells with 1 > s > 0) = 0
    % f(cells with s == 0)    < 0
    f = rldecode(n, nc)-cellNoInCol+1;
    
    % completely filled cells
    s(Gt.columns.cells(f>0)) = 1*(1-sw);
    seff(Gt.columns.cells(f>0)) = 1-sw -sg;
    
    %partially filled cells
    s(Gt.columns.cells(f==0)) = t(n<nc)*(1-sw);
    seff(Gt.columns.cells(f==0)) = t(n<nc)*(1-sw-sg);
    smax = s;
    
    if sg > 0 && any(hmax > h)
        % %hysteresis:
        [n_sr, t_sr] = fillDegree(hmax, Gt);
        
        % remove all cells where hmax - h == 0 and leave only the fraction that is
        % residual co2 in cells that have both residual and free co2
        ix = find(n_sr == n);
        t_sr(ix) = max(t_sr(ix) - t(ix), 0);
        
        ix2 = n_sr - n >= 1;
        f_sr = rldecode(n_sr, nc) - cellNoInCol + 1;
        
        % cells with residual saturation in the whole cell
        ix_tmp = (Gt.columns.cells(f_sr>0 & f<0));
        s(ix_tmp) = sg;
        smax(ix_tmp) = 1-sw;
        
        % cells with residual saturation in bottom part of a cell and free co2 on top
        ix_tmp = Gt.columns.cells(f_sr>0 & f==0);
        currSat = s(ix_tmp);
        s(ix_tmp) = currSat+(1-t(ix2))*sg;
        smax(ix_tmp) = 1-sw;
        
        % cells with possible residual saturation in part of the cell and water in the bottom
        ix_tmp = Gt.columns.cells(f_sr==0);
        currSat = s(ix_tmp);
        s(ix_tmp) = currSat + t_sr(n_sr<nc)*sg;
        smax(ix_tmp) = currSat + t_sr(n_sr<nc) * (1-sw);
    end
end

