function g = oneSlopingFault(dims, drop)
%Make a GRDECL structure for a box grid with a single sloping fault.
%
% SYNOPSIS:
%   grdecl = oneSlopingFault(n, drop)
%
% DESCRIPTION:
%   This function creates a (possibly non-matching) nx-by-ny-by-nz cell
%   corner-point grid containing a single fault.  The resulting geometry
%   discretises the physical domain
%
%      [0,500] x [0,100] x [0,32].
%
% PARAMETERS:
%   n    - Three-element vector, `[nx, ny, nz]`, specifying the number of
%          cells in the `x`, `y`, and `z` coordinate directions
%          respectively.
%
%   drop - Length of fault drop (i.e., the physical length (in metres)
%          by which the two faulted blocks should be offset (along the
%          fault)).  OPTIONAL.  Default value: drop = 0.
%
% RETURNS:
%   grdecl - A `GRDECL` structure suitable for passing to function
%            `processGRDECL`.
%
% NOTE:
%   This example is due to knl.
%
% EXAMPLE:
%   % Create a 90-by-10-by-16 fault grid file with a fault drop of 5 (m)
%   grdecl = oneSlopingFault([90, 10, 16], 5);
%
%   % Create the grid data structure
%   G = computeGeometry(processGRDECL(grdecl));
%
%   % Plot the geometry
%   hg = plotGrid(G, 'FaceAlpha', 0.625, 'EdgeAlpha', 0.1);
%   view(3), grid on, axis tight
%
% SEE ALSO:
%   `processGRDECL`, `writeGRDECL`.

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


% Preprocess input data to make the underlying logical Cartesian grid
if mod(dims(1), 2) ~= 0,
   error('Number of cells in ''x'' direction must be even.');
end
g.cartDims = reshape(dims, 1, []);
[X, Y, Z] = ndgrid(linspace(0, 1, dims(1) + 1), ...
                   linspace(0, 1, dims(2) + 1), ...
                   linspace(0, 1, dims(3) + 1));

X = X + 0.2*(0.5 - abs(X - 0.5)).*(Z - 0.5);
X = X .* 500;
Y = Y .* 100;
Z = Z .*  32;

% Make pilars
lines          = zeros([prod(dims([1, 2]) + 1), 6]);
lines(:,[1,4]) = reshape(X(:,:,[1,end]), [], 2);
lines(:,[2,5]) = reshape(Y(:,:,[1,end]), [], 2);
lines(:,[3,6]) = reshape(Z(:,:,[1,end]), [], 2);

g.COORD = reshape(lines', [], 1);

% Assign z-coordinates
% ind(d) == [1, 2, 2, 3, 3, ..., dims(d), dims(d), dims(d)+1]
ind = @(d) 1 + fix((1 : 2*dims(d)) ./ 2);
z   = Z(ind(1), ind(2), ind(3));

% Add fault
z(end/2+1:end,:,:) = z(end/2+1:end,:,:) + drop;

g.ZCORN = z(:);

% Assign active cells
g.ACTNUM = ones([prod(dims), 1], 'int32');
