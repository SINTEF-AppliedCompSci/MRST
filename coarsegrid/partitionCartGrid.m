function p = partitionCartGrid(cartDims, partDims)
%Partition a Cartesian grid.
%
% SYNOPSIS:
%   partition = partitionCartGrid(cartDims, partDims)
%
% PARAMETERS:
%   cartDims - [nx ny nz] vextor of fine-grid cell dimensions.
%
%   partDims - [cnx cny cnz] vector of coarse-grid block dimensions.
%
% RETURNS:
%   partition - Vector of size [nx*ny*nz 1] with entries equal to
%               coarse block index.
%
% SEE ALSO:
%   `processPartition`.

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


assert(numel(cartDims) == numel(partDims), ...
   'The number of elements in cartDims and partDims differ.');
assert(numel(cartDims) == 2 || numel(cartDims) == 3, ...
   'There must be 2 or 3 positive integers in cartDims,');
assert(all(cartDims > 0) && all(partDims > 0), ...
   'All elements of cartDims and partDims must be positive.');

if numel(cartDims) == 2,
   cartDims = [cartDims(:); 1];
   partDims = [partDims(:); 1];
end

nx  = cartDims(1); ny  = cartDims(2); nz  = cartDims(3);
cnx = partDims(1); cny = partDims(2); cnz = partDims(3);

xC = ceil(reshape(repmat(1:nx, 1    , ny*nz), [], 1) ./ (nx / cnx));
yC = ceil(reshape(repmat(1:ny, nx   ,    nz), [], 1) ./ (ny / cny));
zC = ceil(reshape(repmat(1:nz, nx*ny, 1    ), [], 1) ./ (nz / cnz));

p = xC + cnx*((yC-1) + cny*(zC-1));
