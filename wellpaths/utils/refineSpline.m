function [pts, v] = refineSpline(points, n_refine, interpType)
% Refine a curve to higher resolution using spline interpolation
%
% SYNOPSIS:
%   pts = refineSpline(points, 10, 'spline');
%   pts = refineSpline(points, 5);
%
% DESCRIPTION:
%   Refine a given curve given as a array of points into 
%
% PARAMETERS:
%   points     - A npts x dim array of points giving the curve to be
%                refined. Implicitly assumed to be ordered.
%
%   n_refine   - The refinement factor. If the original entries in points
%                contained n points, the output will have n_refine*n total
%                points.
%
%   interpType - Type of interpolation. Supports the same values as the
%                fourth argument to MATLABs interp1 function. If omitted,
%                it defaults to 'spline'.
%
% RETURNS:
%   pts        - Refined points.
%
%   v          - Parametrization of the new points. Continuous values from
%                1 to npts indicating how far along interpolated values are
%                on the original trajectory.

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

    if nargin == 2
        interpType = 'spline';
    end
    
    n = size(points, 1);
    v0 = (1:n);
  
    dx = (n-1)/(n_refine*n - 1);
    v = (1:dx:n)';
    
    pts = interp1(v0, points, v, interpType);    
end