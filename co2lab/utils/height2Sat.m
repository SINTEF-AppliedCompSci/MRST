function s = height2Sat(sol, g, fluid)
%Convert from height to saturation
%
% SYNOPSIS:
%   s = height2Sat(sol, G, fluid)
%
% PARAMETERS:
%   h - Surface depth.  One scalar value for each column in the top-surface
%       grid.
%
%       Values less than zero are treated as zero while values below the
%       bottom of a column are treated as the column depth.
%
%   G - A top-surface grid as defined by function 'topSurfaceGrid'.
%
%   fluid - a fluid object for example initiated with initVEFluid
%
% RETURNS:
%   s - Saturation - one value for each cell in the underlying 3D model.
%   Corresponds to state.s for the 3D problem.
%
% SEE ALSO:
%   accumulateVertically, integrateVertically

%{
#COPYRIGHT#
%}


if(~isfield(fluid, 'sw'))
   [mu, rho, sr] = fluid.properties();                                     %#ok
   fluid.sr = sr(1);
   fluid.sw = sr(2);
end
% magn = @(v)(sqrt(sum(v.^2,2)));
%n    = g.cells.normals(:,3)./magn(g.cells.normals);
%h = n.*h;
s = zeros(numel(g.columns.cells),1);
% n: number of completely filled cells
% t: fill degree for columns single partially filled cell
[n, t] = fillDegree(sol.h, g); %

% number of cells in each column
nc = diff(g.cells.columnPos);

% compute internal cellNumber in the column for each cell
cellNoInCol = mcolon(ones(g.cells.num,1), diff(g.cells.columnPos))';

% f(cells with s == 1)    > 0
% f(cells with 1 > s > 0) = 0
% f(cells with s == 0)    < 0
f = rldecode(n, nc)-cellNoInCol+1;

% completely filled cells
s(g.columns.cells(f>0)) = 1*(1-fluid.sw);

%partially filled cells
s(g.columns.cells(f==0)) = t(n<nc)*(1-fluid.sw);

if fluid.sr>0 &&any(sol.h_max>sol.h)
   % %hysteresis:
   [n_sr, t_sr] = fillDegree(sol.h_max, g);
   % remove all cells where sol.h_max-sol.h == 0 and leave only the fraction that is
   % residual co2 in cells that have both residual and free co2
   ix = find(n_sr == n);
   t_sr(ix) = max(t_sr(ix)-t(ix),0);

   ix2 = n_sr-n >= 1;
   f_sr = rldecode(n_sr, nc)-cellNoInCol+1;
   %s_copy = s;
   % cells with residual saturation in the whole cell
   s(g.columns.cells(f_sr>0 & f<0)) = fluid.sr;
   % cells with residual saturation in bottom part of a cell and free co2 on top
   currSat = s(g.columns.cells(f_sr>0 &f ==0));
   s(g.columns.cells(f_sr>0 & f==0)) = currSat+(1-t(ix2))*fluid.sr;
   % cells with possible residual saturation in part of the cell and water in the bottom
   currSat = s(g.columns.cells(f_sr==0));
   s(g.columns.cells(f_sr==0)) = currSat + t_sr(n_sr<nc)*fluid.sr;

   %residual = sum(s-s_copy)

end
end

