function grdecl = threeLayers(nx, ny, nz, varargin)
%Construct a corner point discretization of a three-layered structure.
%
% SYNOPSIS:
%   grdecl = threeLayers(nx, ny, nz)
%
% DESCRIPTION:
%   This function creates a corner point discretization of a mapped
%   three-layered structure from the physical domain
%
%      [-1,1] x [0,1] x [0,1].
%
%   The horizontal regions correspond to dividing planes at z=0.25 and
%   z=0.75.  Moreover, when x>=0, the z-coordinate of each vertex is mapped
%   according to the rule
%
%      z = z - (z + 1)*y/2,  y \in [0, 1]
%
% PARAMETERS:
%
%   nx - One- or two-element vector defining the number of intervals along
%        the X axis.  Specifically, if nx=[n1, n2], then the interval
%        [-1,0] will be divided into n1 sub-intervals and the interval
%        [0,1] will be divided into n2 sub-intervals.  A scalar `nx`
%        parameter is treated as if the caller specified nx = [nx, nx].
%
%   ny - Number of sub-intervals along the Y axis.  Single positive integer.
%
%   nz - One or three-element vector defining the number of sub-intervals
%        along the Z axis.  Specifically, if `nz=[n1, n2, n3]`, then the
%        interval `[0,0.25]` will be divided into `n1` sub-intervals, the
%        interval `[0.25,0.75]` will be divided into `n2` sub-intervals,
%        and the interval `[0.75,1]` will be divided into `n3`
%        sub-intervals. A scalar `nz` parameter is treated as if the caller
%        specified nz = [nz, nz, nz].
%
% EXAMPLE:
%   grdecl = threeLayers(10, 10, [5, 10, 5]);
%   G = processGRDECL(grdecl)
%   h = plotGrid(G, 'EdgeAlpha', 0.25, 'FaceAlpha', 0.375);
%   set(get(h, 'Parent'), 'ZDir', 'normal')
%   view(148.5, 32)
%
% RETURNS:
%   grdecl - A GRDECL structure suitable for further processing by
%            function 'processGRDECL' or for permanent storage on disk by
%            function 'writeGRDECL'.
%
% SEE ALSO:
%   `processGRDECL`, `writeGRDECL`, `simpleGrdecl`, `makeModel3`.

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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


   if numel(nx) == 1, nx = [nx, nx]    ; else assert (numel(nx) == 2), end
                                              assert (numel(ny) == 1);
   if numel(nz) == 1, nz = [nz, nz, nz]; else assert (numel(nz) == 3), end

   % Define coordinates.
   x = [linspace(-1   , 0   , nx(1) + 1), ...
        linspace( 0   , 1   , nx(2) + 1)];      x = unique(x);
   y =  linspace( 0   , 1   , ny    + 1);
   z = [linspace( 0   , 0.25, nz(1) + 1), ...
        linspace( 0.25, 0.75, nz(2) + 1), ...
        linspace( 0.75, 1   , nz(3) + 1)];      z = unique(z);

   [X, Y, Z] = ndgrid(x, y, z);

   % DIMENS (cartDims).
   grdecl.cartDims = [numel(x), numel(y), numel(z)] - 1;

   % COORDS.
   nlines = prod(grdecl.cartDims(1:2) + 1);
   lines  = zeros([nlines, 6]);
   lines(:, [1, 4]) = reshape(X(:,:,[1, end]), [], 2);
   lines(:, [2, 5]) = reshape(Y(:,:,[1, end]), [], 2);
   lines(:, [3, 6]) = reshape(Z(:,:,[1, end]), [], 2);
   grdecl.COORD = reshape(lines.', [], 1);

   % ZCORN.
   % ind(d) == [1, 2, 2, 3, 3, ..., dims(d), dims(d), dims(d)+1]
   ind = @(d) 1 + fix((1 : 2*grdecl.cartDims(d)) ./ 2);
   y   = Y(ind(1), ind(2), ind(3));
   z   = Z(ind(1), ind(2), ind(3));
   i   = rldecode([false; true], grdecl.cartDims([1; 1]));
   z(i,:,:) = z(i,:,:) - (z(i,:,:) + 1).*y(i,:,:)./2;

   grdecl.ZCORN  = reshape(z, [], 1);
   grdecl.ACTNUM = ones([prod(grdecl.cartDims), 1], 'int32');

   % PERMX.
   cnz   = cumsum(nz);
   permx = zeros(grdecl.cartDims);
   permx(:,:, 0      + 1 : cnz(1)) = 100;
   permx(:,:, cnz(1) + 1 : cnz(2)) =   1;
   permx(:,:, cnz(2) + 1 : cnz(3)) = 100;

   grdecl.PERMX = reshape(permx, [], 1);
end
