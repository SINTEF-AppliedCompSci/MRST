function [A,n] = rlencode(A, dim)
%Compute run length encoding of array A along dimension dim.
%
% SYNOPSIS:
%   [A,n] = rlencode(A)
%   [A,n] = rlencode(A, dim)
%
% PARAMETERS:
%   A         - Array
%   dim       - dimension of `A` where run length encoding is done.
%               `dim > 0`, `dim <= ndims(A)`.
%               OPTIONAL.  Default value: `dim=1`.
%
% RETURNS:
%   A         - Compressed `A` where repeated layers are removed.
%   n         - repetition count of repeated layers in original `A`.
%
% EXAMPLE:
%
%   A = [1,2,3,4;1,2,3,4;3,4,5,6;3,3,3,3;3,3,4,5;3,3,4,5]
%   [A,n] = rlencode(A,1)
%
%
% SEE ALSO:
%   `rldecode`

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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


if nargin < 2
  dim = 1;
end

if isempty(A), n = 0; return; end

% Take dimension we compress along to be first dimension,
% i.e., swap dimensions 1 and dim.

d      = 1:ndims(A);
d(1)   = dim;
d(dim) = 1;
B      = permute(A,d);

% Pick out (1:end,:,...) and (2:end,:,...) in
% a multidimensional way


% Find positions where layers differ
nanB = isnan(B);
%i    = [find(any(B(1:end-1,:) ~= B(2:end,:)), 2); size(B,1)];
i    = [find(any( xor((B(1:end-1,:) ~= B(2:end,:)) , ...
                       nanB(1:end-1,:)==1 & nanB(2:end,:)==1), 2)); size(B,1)];
% compare differences in position to find run length.
n = diff([0;i]);

% swap dimensions 1 and dim.
sz = size(B);sz(1)=numel(i);
A  = permute(reshape(B(i,:), sz),d);
end
