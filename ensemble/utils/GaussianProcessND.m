function res = GaussianProcessND(N, corr_fun)
%
% Simulate an d-dimensional Gaussian process for a given covariance function
%
% SYNOPSIS:
%   function res = GaussianProcessND(N, corr_fun)
%
% DESCRIPTION:
% 
%  The result is computed on the [0,1]^d box, with a number of samples specified
%  by the n-dimensional vector 'N'. The covariance function function 'g' is a
%  nonnegative function from [0,1] to [0,1], converting a distance between two
%  points on the [0,1]^d box to a covariance value.  In general, this
%  should be a decreasing function such that g(0) = 1.  Well-known covariance
%  functions are the exponential covariance function, the squared exponential
%  covariance function and the spherical covariance function.
% 
% PARAMETERS:
%   N        - d-dimensional vector giving the number of samples along each direction
%   corr_fun - covariance function (nonnegative from [0,1] to [0,1]).
%
% RETURNS:
%   res - a realization of the d-dimensional Gaussian process.
%
% EXAMPLES:
%
% 1) Simulate a 2D field with 100 samples in each direction, using an exponential
%    covariance function (yields an irregular field):
% 
%    res = GaussianProcessND([100, 100], @(xy) exp(-sqrt(sum(xy.^2,2))/0.3));
%  
% 2) Simulate a 2D field with 100 samples in each direction, using a squared
%    exponential covariance function (yields a smooth field):
% 
%    res = GaussianProcessND([100, 100], @(xy) exp(-sum(xy.^2,2)*5e2))  
%  
% SEE ALSO:
% GaussianProcess1D - the one-dimensional version of this function

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

  % Determine size and eigenvalues of the block-circulant covariance matrix C
  lambda = determine_circular_covariance_matrix_eigenvalues(N, corr_fun);
  
   % Compute W = \lambda^{-1/2} Q^* Z 
   % (See p. 413 of Wood, A. and Chan, G. (1994)) (*) 
   W = simulate_half_transformed_gaussian_variable(lambda);
   
   % Apply the full transformation on W (i.e. the product Q W) to obtain a real
   % gaussian vector governed by the provided correlation function.  In
   % practice, we obtain this using a fast fourier transform, and keeping
   % only the N first values (the rest are only present to satisfy the
   % construction of the circulant matrix).
   res = real(fftn(W));  % real gets rid of remaining imaginary numerical noise
   
   clip = cell(numel(N), 1);
   for d = 1:numel(N)
      clip{d} = 1:N(d);
   end
   res = res(clip{:});

end

% (*) "Simulation of Stationary Gaussian Processes in [0, 1]^d", Andrew Wood
% and Grace Chan, Journal of Computational and Graphical Statistics, Vol.3
% Number 4 (1994), pages 409-432

% ----------------------------------------------------------------------------
function c = compute_matrix_coefficients(N, g, corr_fun)

   m = 2.^g; 
   ixs = arrayfun(@(ix) (-ix+1):(ix-1), m, 'UniformOutput', false);
   hix = cell(ixs);
   [hix{:}] = ndgrid(ixs{:});
   hix = cellfun(@(m) m(:), hix, 'UniformOutput', false);
   h = [hix{:}];
   
   arg = h; % argument to correlation function, equals h_tilde / N according
            % to Wood and Chan (1994).

   for d = 1:numel(N)
      high_ix = arg(:,d) > m(d)/2;
      low_ix = -arg(:,d) > m(d)/2;
      
      arg(high_ix, d) = arg(high_ix, d) - m(d);
      arg(low_ix , d) = arg(low_ix, d)  + m(d);
      
      arg(:,d) = arg(:,d) ./ N(d);
   end

   c = corr_fun(arg);
   
   % preventing potential problem for 'uneven' correlation functions
   % (cf. Wood and Chan (1994) page 415)
   for d = 1:numel(N)
      c(abs(h(:,d)) == m(d)/2) = 0;
   end

   % the construction of c is now finished, and it should now have a block
   % circulant structure.  We keep only the part we care about (for which we
   % are later going to apply the fast fourier transform).
   c = reshape(c, [2.^(g+1)-1, 1]);
   clip = cell(numel(N), 1);
   for d = 1:numel(N)
      clip{d} = ceil(size(c, d)/2):size(c, d);
   end
   c = c(clip{:});
end

% ----------------------------------------------------------------------------
function lambda = determine_circular_covariance_matrix_eigenvalues(N, corr_fun)

   D = numel(N); % number of dimensions
   gmin = ceil(log2(2 * (N-1)));
   EPS = 1e-10;

   % increasing g until all eigenvalues of C are positive (or until the
   % matrix has grown so big that we give up)
   g_inc = 10;  % we allow matrix to double in size this number of times
   g = gmin;
   found = false;
   for i = 1:g_inc
      c = compute_matrix_coefficients(N, g, corr_fun);
      
      lambda = real(fftn(c)); % components should be real, but we apply
                              % 'real()' to remove any imaginary numerical
                              % residual 
      if min(lambda(:)) >= -EPS
         lambda(lambda < 0) = 0; % numerical noise
         found = true;
         break;
      end
      g(mod(i-1, D)+1) = g(mod(i-1,D)+1) + 1;
   end
   
   if ~found
      error('Could not find a positive semi-definite circulant matrix.');
   end
   
   % attempt to decrease g while keeping eigenvalues of C positive
   while(true)
      reducible_dimension = 0; 
      for dim = 1:D
         % try to reduce size along dimension 'dim'
         if g(dim) > gmin(dim)
            g_tmp = g;
            g_tmp(dim) = g(dim) - 1;
            c_tmp = compute_matrix_coefficients(N, g_tmp, corr_fun);
            lambda_tmp = real(fftn(c_tmp));
            if min(lambda_tmp(:)) >= -EPS
               reducible_dimension = dim;
               lambda = lambda_tmp; % keep the eigenvals for the smaller matrix
               break;
            end
         end
      end
      
      if reducible_dimension > 0 
         g(reducible_dimension) = g(reducible_dimension) - 1;
      else
         % no dimension could be further reduced.  Give up further reduction
         break; 
      end      
   end
   
end

% ----------------------------------------------------------------------------
function W = simulate_half_transformed_gaussian_variable(lambda)
   
   ixs = arrayfun(@(ix) 1:ix, size(lambda), 'UniformOutput', false);
   hix = cell(ixs);
   [hix{:}] = ndgrid(ixs{:});
   hix = cellfun(@(x) x(:), hix, 'UniformOutput', false);
   h = [hix{:}];

   % for each entry in h, identify the entry for which its stochastic
   % variable is correlated.
   for d = 1:size(h,2)
      col = h(:,d);
      modif_ix = logical((col ~= 1) & (col ~= size(lambda, d)/2+1));
      h(modif_ix, d) = size(lambda, d) + 2 - h(modif_ix, d);
   end

   % for each row of 'h', 'I' will give us the row corresponding to the
   % single stochastic variable correlated with this one.
   [~, I] = sortrows(h, numel(size(lambda)):-1:1);

   W = zeros(numel(I), 1);
   m = numel(lambda);

   I0 = (1:numel(I))';
   same_ix = (I0 == I); num_same = sum(same_ix);
   lower_ix = (I0 < I); num_lower = sum(lower_ix);

   % rows pointing to themselves are assigned independent random
   % variables
   W(same_ix) = sqrt(lambda(same_ix)/m) .* randn(num_same, 1);
   
   % Other row-pairs are assigned correlated (1/2) random variables
   U = randn(num_lower, 1); V = randn(num_lower, 1);
   Za = U + 1i * V;  Zb = U - 1i * V;
   fac = sqrt(lambda(lower_ix)/(2*m));

   W(lower_ix) =  fac .* Za;
   W(I(lower_ix)) = fac .* Zb;
   
   W = reshape(W, size(lambda));
end
