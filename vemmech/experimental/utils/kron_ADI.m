function res = kron_ADI(u, v)%
%
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
%

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
