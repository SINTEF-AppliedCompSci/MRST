function res = kron_ADI(u, v)
% A version of kron that supports ADI vectors
%
% SYNOPSIS:
%   function res = kron_ADI(u, v)
%
% DESCRIPTION:
%
% PARAMETERS:
%   u - vector or ADI vector
%   v - vector or ADI vector
%
% RETURNS:
%   res - kron product

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

% only supports vectors, not matrices
assert(isa(u, 'ADI') || isvector(u));
assert(isa(v, 'ADI') || isvector(v));

numel_u = numel(u);
if isa(u, 'ADI')
   numel_u = numel(u.val);
end

numel_v = numel(v);
if isa(v, 'ADI')
   numel_v = numel(v.val);
end

ix = (1:numel_u * numel_v)';
iy = repelem((1:numel_u)', numel_v, 1);


mat = sparse(ix, iy, 1, numel_u * numel_v, numel_u);

res = (mat * u) .* repmat(v, numel_u, 1);
res = full(res);
