function res = GaussianProcess1D(N, corr_fun) 
%
% Simulate a one-dimensional Gaussian process for a given covariance
% function on the [0, 1] interval.
% 
% SYNOPSIS:
%   function res = GaussianProcess1D(N, corr_fun)
%
% DESCRIPTION:
% 
%  The result is computed on the [0,1] interval, with a number of samples
%  specified by the argument 'N'. The covariance function function 'g' is a
%  nonnegative function from [0,1] to [0,1], converting a distance between two
%  points on the [0, 1] interval to a correllation vlaue.  In general, this
%  should be a decreasing function such that g(0) = 1.  Well-known
%  covariance  functions are the exponential covariance function, the squared
%  exponential covariance function and the spherical covariance function.
% 
% PARAMETERS:
%   N        - Number of samples on the [0,1] interval.
%   corr_fun - variance function (nonnegative from [0,1] to [0,1]).
%
% RETURNS:
%   res - a realization of the Gaussian process.
%
% EXAMPLE:
%
% 1) Simulate a Gaussian process using 1000 samples and an exponential
%    covariance function (yields an irregular field):
% 
%    res = GaussianProcess1D(1000, @(x) exp(-sqrt(x.^2)/0.3));
%  
% 2) Simulate a Gaussian process using 1000 samples and a squared
%    exponential covariance function (yields a smooth field):
% 
%    res = GaussianProcess1D(1000, @(x) exp(-x.^2*5e2))  
%
% SEE ALSO:
% GaussianProcessND - multivariate version of this function

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

   % Determine size and eigenvalues of the circulant covariance matrix C.  
   lambda = determine_circulant_covariance_matrix_eigenvalues(N, corr_fun);
   
   % Compute W = \lambda^{-1/2} Q^* Z 
   % (See p. 413 of Wood, A. and Chan, G. (1994)) (*) 
   W = simulate_half_transformed_gaussian_variable(lambda);
   
   % Apply the full transformation on W (i.e. the product Q W) to obtain a real
   % gaussian vector governed by the provided correlation function.  In
   % practice, we obtain this using a fast fourier transform, and keeping
   % only the N first values (the rest are only present to satisfy the
   % construction of the circulant matrix).
   res = fft(W);
   res = real(res(1:N)); % we remove any numerical residual of imaginary components
end

% (*) "Simulation of Stationary Gaussian Processes in [0, 1]^d", Andrew Wood
% and Grace Chan, Journal of Computational and Graphical Statistics, Vol.3
% Number 4 (1994), pages 409-432

% ----------------------------------------------------------------------------
function lambda = determine_circulant_covariance_matrix_eigenvalues(N, corr_fun)
   
   % determine 'g' such that 2^g >= 2 * (N-1) and that the circulant
   % covariance matrix C(2^g, 2^g) is positive semi-definite
   gmin = ceil(log2(2 * (N-1)));
   gmax = gmin + 10;
   EPS = 1e-10;
   
   for g = gmin:gmax
      vec = zeros(1, 2^g);
      half_ix = 2^(g-1);
      vec(1:half_ix+1) = corr_fun((0:half_ix)'./N);
      vec(half_ix+2:end) = fliplr(vec(2:half_ix));
      
      lambda = fft(vec);
      
      % Due to symmetry, 'lambda' should contain real values.  However, for
      % numerical reasons there may be residual imaginary components.  We
      % wrap 'lambda' in a call to 'real' below when determining whether all
      % eigenvalues are positive.
      if min(real(lambda)) >= -EPS
         % all eigenvalues are positive or zero (disregarding numerical noise).  We can
         % use this circulant matrix to embed our covariance matrix.
         lambda(lambda < 0) = 0; % numerical noise
         return;
      end
   end

   % If we got here, we did not succeed in finding a positive semi-definite
   % circulant matrix
   error('Did not succeed in finding a positive semi-definite circulant matrix.');
   
end


% ----------------------------------------------------------------------------
function W = simulate_half_transformed_gaussian_variable(lambda)
   
   W = zeros(size(lambda));
   m = numel(lambda);
   mid = m/2; % m should be a power of two, so 'mid' should be a whole number
   lambda_sqrt = sqrt(real(lambda)); % use 'real' to get rid of possible
                                     % residual imaginary component

   % Fill in the two completely uncorrelated random variables of the
   % transformed system
   W(1) = lambda_sqrt(1) * randn();
   W(mid+1) = lambda_sqrt(mid+1) * randn();
   
   % Fill in the pairwise correlated random variables
   U = randn(1, mid-1);
   V = randn(1, mid-1);
   
   W(2:mid) = 1/sqrt(2) * lambda_sqrt(2:mid) .* (U + 1i * V);
   
   tmp = 1/sqrt(2) * lambda_sqrt(2:mid) .* (U - 1i * V);
   W(mid+2:end) = fliplr(tmp);
   
   % Scale W correctly
   W = W / sqrt(m);
   
end