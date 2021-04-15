function [K, L] = logNormLayers(N, varargin)
%Compute realization of lognormal, isotropic permeability field.
%
% SYNOPSIS:
%   K     = logNormLayers(N)
%   [K,L] = logNormLayers(N, M)
%   [K,L] = logNormLayers(N, M, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function creates a nx-by-ny-by-nz scalar permeability field
%   consisting of `M` layers. If `M` is a vector, then the number of layers
%   equals `numel(M)` and `M` gives the desired mean for each layer. If `L`
%   is given, we assume that `M` gives the desired mean values of each
%   layer and then `numel(M)=numel(L)-1`.
%
%   The gaussian field within each layer is generated based on the
%   formulas::
%     k = smooth3(randn(..) - a*randn(1), 'gaussian', sz, std);
%     k = exp( b + sigma*k);
%   where `a`, `b`, `sigma`, `sz`, and `std` are optional parameters.
%
% PARAMETERS:
%   N       - Three-element vector, `[nx, ny, nz]`, specifying the number of
%             discrete values in the `x`, `y`, and `z` coordinate
%             directions respectively.
%
%   M      -  The number of layers if a scalar and the desired means of the
%             layers if a vector.
%
%
% OPTIONAL PARAMETERS:
%
%   'Verbose' - Whether or not to emit progress reports during the
%               computation. 
%               Logical.  Default value: `mrstVerbose()`.
%
%   'indices' - List `L` of k-indices for the generated layers, given such
%               that layer <i> is represented by k-indices `L(i):L(i+1)-1`
%               in the grid. The number of layers must agree with the
%               specification of means in the `M` parameter.
%               Default value: `[]`
%
%    'sigma   - Scalar. Default value: `2`
%
%    'a'      - Scalar. Default value: `0.6`
%
%    'b'      - Scalar. Default value: `2`
%
%    'sz'     - Gaussian mask. Default value: `[9,3,3]`
%
%    'std'    - Standard deviation. Default value: `4.5`
%
% RETURNS:
%
%   K - The scalar nx-by-ny-by-nz permeability field in units mD
%
%   L - list of indices for the layers in `K`, so that layer `i` has
%       k-indices=`L(i):L(i+1)-1`

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


opt = struct('Verbose', mrstVerbose, ...
             'indices', [], ...
             'sigma',   2, ...
             'a',       0.6, ...
             'b',       2, ...
             'sz',      [9, 3, 3], ...
             'std',     4.5);

if all(size(N) == [1,2]) && mod(nargin, 2)
   N = [N, 1];
   M = 1;
else
   if nargin==1
      M = 1;
   else
      M = varargin{1};
   end
   varargin = varargin(2:end);
end
opt = merge_options(opt, varargin{:});
rand_layers = @(m,n) unique([1, fix(1 + rand([1,m-1])*(n-1)), n+1]);

% Process input parameters
if ~isempty(opt.indices)
   L = opt.indices;
   if (max(L)~=N(3)+1) || any(diff(L)<=0) || any(L<0)
      error(msgid('NumLayers:Iconsistent'), ...
         'Inconsistent values. L must increase in [1,%d]', N(3)+1);
   elseif numel(L)>N(3)+1
      error(msgid('NumLayers:Excessive'), ...
           ['Number of requested permeability layers (%d) ', ...
            'exceeds number of model layers (%d).'], L(end), N(3));
   elseif numel(L)~=numel(M)+1
      error(msgid('NumLayers:Inconsistent'),...
         'Incorrect input: numel(M) should be equal numel(L)-1');
   end
   lmean = M(:);
elseif (nargin < 2) || isempty(M)
   L = [1:N(3), N(3)+1];  lmean = [];
elseif numel(M) == 1
   % i.e., only the number of layers is specified
   if M > N(3)
      error(msgid('NumLayers:Excessive'), ...
           ['Number of requested permeability layers (%d) ', ...
            'exceeds number of model layers (%d).'], M, N(3));
   end
   if M < 0
      L = [1, N(3)+1];
   elseif M == N(3)
      L = 1 : N(3) + 1;
   else
      L = [];
      while numel(L) < M + 1
         L = rand_layers(M, N(3));
      end
   end
   lmean = [];
else
   % i.e., a mean is specified per layer
   if numel(M) > N(3)
      error(msgid('NumLayers:Excessive'), ...
           ['Number of requested permeability layers (%d) ', ...
            'exceeds number of model layers (%d).'], numel(M), N(3));
   end
   if numel(M) == N(3)
      L = 1 : N(3) + 1;
   else
      L = [];
      while numel(L) < numel(M) + 1
         L = rand_layers(numel(M), N(3));
      end
   end
   lmean = M(:);
end

% check input
if numel(opt.a)~=1
   error(msgid('a:wrong'), 'Optimal parameter <a> not scalar');
end
if numel(opt.b)~=1
   error(msgid('b:wrong'), 'Optimal parameter <b> not scalar');
end
if numel(opt.sigma)~=1
   error(msgid('std:wrong'), 'Optimal parameter <sigma> not scalar');
end

% Generate layers
if opt.Verbose
   fprintf('\nGenerating lognormal, layered permeability\n  Layers: [ ');
end
for i = 1 : length(L)-1
    n = L(i+1) - L(i);
    k = smooth3(randn([N(1:2), n+2]) - opt.a*randn(1), ...
                'gaussian', opt.sz, opt.std);
    k = exp(opt.b + opt.sigma*k);
    k = k(:,:,2:n+1);
    if numel(lmean) > 0
       K(:,:,L(i):L(i+1)-1) = lmean(i)*k / mean(k(:));
    else
       K(:,:,L(i):L(i+1)-1) = k;
    end
    if opt.Verbose, fprintf('%d:%d ', L(i), L(i+1)-1); end
end
K = K(:);

if opt.Verbose
   fprintf(']\n');
   fprintf('  min: %g, max: %g [mD], ratio: %g\n\n', ...
           min(K), max(K), max(K) / min(K));
end
