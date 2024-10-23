function s = height2Sat(sol, Gt, fluid)
%Convert from height to saturation
%
% SYNOPSIS:
%   s = height2Sat(sol, Gt, fluid)
%
% PARAMETERS:
%   h - CO2 plume thickness.  One scalar value for each column in the
%       top-surface grid.
%
%       Values less than zero are treated as zero while values below the
%       bottom of a column are treated as the column depth.
%
%   Gt - A top-surface grid as defined by function 'topSurfaceGrid'.
%
%   fluid - a fluid object for example initiated with initVEFluid
%
% RETURNS:
%   s - Saturation - one value for each cell in the underlying 3D model.
%   Corresponds to state.s for the 3D problem.
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

deprecatedMessage('height2Sat', 'height2finescaleSat');

if isfield(fluid, 'res_gas')
   % ADI-type fluid
   sr = fluid.res_gas;
   sw = fluid.res_water;
elseif(~isfield(fluid, 'sw'))
   [mu, rho, res] = fluid.properties();                                     %#ok
   sr = res(1);
   sw = res(2);
end

% magn = @(v)(sqrt(sum(v.^2,2)));
%n    = Gt.cells.normals(:,3)./magn(Gt.cells.normals);
%h = n.*h;
s = zeros(numel(Gt.columns.cells),1);
% n: number of completely filled cells
% t: fill degree for columns single partially filled cell
[n, t] = fillDegree(sol.h, Gt); %

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

%partially filled cells
s(Gt.columns.cells(f==0)) = t(n<nc)*(1-sw);

if sr>0 &&any(sol.h_max>sol.h)
   % %hysteresis:
   [n_sr, t_sr] = fillDegree(sol.h_max, Gt);
   % remove all cells where sol.h_max-sol.h == 0 and leave only the fraction that is
   % residual co2 in cells that have both residual and free co2
   ix = find(n_sr == n);
   t_sr(ix) = max(t_sr(ix)-t(ix),0);

   ix2 = n_sr-n >= 1;
   f_sr = rldecode(n_sr, nc)-cellNoInCol+1;
   %s_copy = s;
   % cells with residual saturation in the whole cell
   s(Gt.columns.cells(f_sr>0 & f<0)) = sr;
   % cells with residual saturation in bottom part of a cell and free co2 on top
   currSat = s(Gt.columns.cells(f_sr>0 &f ==0));
   s(Gt.columns.cells(f_sr>0 & f==0)) = currSat+(1-t(ix2))*sr;
   % cells with possible residual saturation in part of the cell and water in the bottom
   currSat = s(Gt.columns.cells(f_sr==0));
   s(Gt.columns.cells(f_sr==0)) = currSat + t_sr(n_sr<nc)*sr;

   %residual = sum(s-s_copy)

end
end

% ----------------------------------------------------------------------------
function [n, t] = fillDegree(h, G)
%Compute degree of fill (saturation).
%
% SYNOPSIS:
%   [n, t] = fillDegree(h, G)
%
% PARAMETERS:
%   h - CO2 plume thickness.  One scalar value for each column in the
%       top-surface   grid.  As a special case a single, scalar value is
%       repeated for each column in the top-surface grid.
%
%       Values less than zero are treated as zero while values below the
%       bottom of a column are treated as the column depth.
%
%   G - A top-surface grid as defined by function 'topSurfaceGrid'.
%
% RETURNS:
%   n - Number of completely filled cells.  One non-negative integer for
%       each grid column.  If the top-most cell in column 'i' is partially
%       filled, then n(i)==0.
%
%   t - Fill degree of column's single partially filled cell.  One scalar
%       value between zero and one (inclusive, i.e., t \in [0,1]) for each
%       column.  If a column is completely filled (i.e., if n(i) equals the
%       number of cells in column 'i'), then t(i) is (arbitrarily) set to
%       zero.
%
% NOTE:
%   Assume column 'i' consists of N(i) cells.  Then, implicitly, the
%   saturation is s==1 in the n(i) first cells of this column.
%   Furthermore, s=t(i) in cell n(i)+1, and s==0 in the remaining
%   N(i) - (n(i) + 1) cells below the single, partially filled cell.
    
    if numel(h) == 1,
        % Scalar h => *same* height in all columns.
        h = h(ones([G.cells.num, 1]));
    end
    
    % Exclude non-physical h<0 case (and guarantee col-vec shape).
    h = reshape(max(h, 0), [], 1);
    
    % Compute look-up table.
    %   - f(c,1): Dual purpose.
    %             i)   f(c,1) < 0 in completely filled cells, 'c'.
    %             ii) -f(c-1,1)/(f(c,1)-f(c-1,1)) == fill degree in cell 'c'
    %                  (and *only* in cell 'c') if f(c-1,1) * f(c,1) < 0.
    %
    %   - f(c,2): Column index of cell 'c'.
    %
    nc     = reshape(diff(G.cells.columnPos), [], 1);  assert (all(nc > 0));
    f      = rldecode([h, (1:G.cells.num)'], nc);
    f(:,1) = G.columns.z - f(:,1);
    
    % Count completely filled cells in each column.
    % n==0 if first cell is partially filled.
    %
    n = accumarray(f(:,2), f(:,1) <= 0, size(h));
    t = zeros(size(n));
    
    % Compute position/row index into look-up table (f) of each column's
    % partially filled cell (out of bounds if n==nc).
    p = G.cells.columnPos(1 : end-1) + n;
    
    % First cell partially filled is a special case.
    % Use the column's h-value directly to compute fill degree.
    %
    first    = n == 0;
    t(first) = h(first) ./ G.columns.z(p(first));
    
    % General case (partially filled column).
    % Formula ii) defines fill degree.
    % Completely filled columns (n==nc) excluded (excess fill degree zero).
    %
    j    = (n < nc) & ~first;
    pj   = p(j);
    t(j) = -f(pj-1, 1) ./ (f(pj, 1) - f(pj-1, 1));
    
    assert (~any((t < 0) | (t > 1)));
end
