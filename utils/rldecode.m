function A = rldecode(A, n, dim)
%Decompress run length encoding of array `A` along dimension `dim`.
%
% SYNOPSIS:
%   B = rldecode(A, n, dim)
%   B = rldecode(A, n) % dim assumed to be 1
%
% PARAMETERS:
%   A         - encoded array
%   n         - repetition of each layer along dimension `dim`. `n` can be
%               either a scalar or one repetition number for each layer.
%   dim       - dimension of `A` where run length encoding is done.
%               dim > 0.
%
% RETURNS:
%   B         - uncompressed `A`
%
% EXAMPLE:
%   % 1) Numerical example:
%   A = [1,2,3,4;1,2,3,4;3,4,5,6;3,3,3,3;3,3,4,5;3,3,4,5]
%   [B,n] = rlencode(A,1)
%   C = rldecode(B,n,1)
%
%   % 2) Retrive 'first' column of G.cells.faces (see grid_structure):
%   G = cartGrid([10, 10, 2]);
%   cellNo  = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
%   disp(['CellFace nr. 10 belongs to cell nr: ', num2str(cellNo(10))]);
%
% SEE ALSO:
%   `rlencode`

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


if nargin < 3
  dim = 1;
end

assert(dim > 0, 'Third argument DIM must be positive');

if numel(n) == 1,
   n = repmat(n, [size(A, dim), 1]);
end

assert (all( n(:)>=0 ), 'All repeat counts should be nonnegative.');
if nargin < 3,
   assert (numel(n) == size(A, dim), ...
   sprintf(['There should be a repeat count for each value along dimension dim.\n',...
    'The default value of dim is 1. Did you forget to specify dim?']));
else
   assert (numel(n) == size(A, dim), ...
   'There should be a repeat count for each value along dimension dim.');
end

% Take dimension we compress along to be first dimension,
% i.e., swap dimensions 1 and dim.
d      = 1:max(dim, ndims(A));
d([1, dim])   = [dim, 1];
B      = permute(A,d);

r      = n(:)~=0;
B      = reshape(B(r, :), sum(r), []);



% Insert repeated layers and permute back
i      = cumsum([1; double(reshape(n(r), [], 1))]);
j      = zeros(i(end)-1,1);
j(i(1:end-1)) = 1;

szA    = [size(A), ones(1, dim-ndims(A))];
A      = permute(reshape(B(cumsum(j),:), [sum(n(:)), szA(d(2:end))]), d);
end
