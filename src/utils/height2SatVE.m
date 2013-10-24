function s = height2SatVE(sol, g, fluid)
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
s=sol.h*(1-fluid.sw)+(sol.h-sol.h_max)*fluid.sr;
s=s./g.cells.H;
end

