function g = makeModel3(dims, varargin)
%Build a synthetic geometry with two faults.
%
% SYNOPSIS:
%   grdecl = makeModel3(dims)
%   grdecl = makeModel3(dims, physDims)
%
% PARAMETERS:
%  dims     - A 3-vector `[nx, ny, nz]` giving the number of cells in each of
%             the three coordinate directions.
%
%  physdims - A 3-vector `[hx, hy, hz]` giving the physical extent (in units
%             of meters) of the original sandbox before the transformations
%             that emulate geological activity. OPTIONAL.  Default value:
%             `[1000, 600, 15]`,  meaning the geometry discretises a sandbox
%             of the following size in physical domain:
%
%                    `[0,1000]` by `[0,600]` by `[0,15]`
%
% RETURNS:
%   grdecl - A GRDECL structure suitable for passing to function
%            `processGRDECL`.
%
% EXAMPLE:
%   plotGrid(processGRDECL(makeModel3([100, 60, 15])));
%
% SEE ALSO:
%   `processGRDECL`

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


physDims = [1000, 600, 15];
if (nargin > 1) && isnumeric(varargin{1}) && (numel(varargin{1}) == 3),
   physDims = varargin{1};
end
g.cartDims = reshape(dims, 1, []);
[X, Y, Z]  = ndgrid(linspace(0, physDims(1), dims(1) + 1), ...
                    linspace(0, physDims(2), dims(2) + 1), ...
                    linspace(0, physDims(3), dims(3) + 1));

% Make perturbation in z-direction
a   = rand ([dims(3) + 1, 1]);
s   = randn([dims(3) + 1, 2]);
fun = @(x, y, a, s) -5 *sin(pi * (x - 1/2)) - ...
                     3 *sin(pi * (y + 2*x)) + ...
                    a/2*sin(2*pi*x + s(1)) .* sin(4*pi*y + s(2));

for k = 1 : dims(3) + 1,
   xi       = X(:,:,k) ./ physDims(1);
   eta      = Y(:,:,k) ./ physDims(2);
   Z(:,:,k) = Z(:,:,k) -  fun(xi, eta, a(k), s(k,:));
end


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

% Add faults
h  = physDims ./ dims;
i1 = round((0.6 * physDims(1)) / h(1));   assert (i1 > 0);
j1 = round((0.5 * physDims(2)) / h(2));   assert (j1 > 0);
i2 = round((0.5 * physDims(1)) / h(1));   assert (i2 > 0);
j2 = round((1/3 * physDims(2)) / h(2));   assert (j2 > 0);

z(1:2*i1,     1:2*j1,     :) = z(1:2*i1,     1:2*j1,     :) - 2;
z(2*i2+1:end, 2*j2+1:end, :) = z(2*i2+1:end, 2*j2+1:end, :) + 3;

% Pinch layers with negative thickness
z = cumsum(cat(3, z(:,:,1), max(0, diff(z, 1, 3))), 3);

g.ZCORN = z(:);

% Assign active cells
actnum = ones(dims);
[x, y] = ndgrid(linspace(0, 1, dims(1)), ...
                linspace(0, 1, dims(2)), ...
                linspace(0, 1, dims(3)));

actnum((x - 0.4).^2 +  y     .^2 > 0.95) = 0;
actnum((x - 0.6).^2 + (y - 1).^2 > 0.9 ) = 0;

g.ACTNUM = int32(actnum(:));
