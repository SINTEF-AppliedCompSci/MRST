function [params, history] = mle_prob_param(samples, distr_fun, distr_fun_grad, init_p, pmin, pmax) 
%
% Maximum likelihood estimation of probability function parameters based on a
% number of measured samples.
% 
% SYNOPSIS:
%   function [params, history] = mle_prob_param(samples, distr_fun, distr_fun_grad, init_p, pmin, pmax)
%
% DESCRIPTION:
% Based on a number of samples ('samples'), which are supposedly drawn from a
% known probability distribution with unknown parameter(s), estimate the
% parameters that would yield the maximum probability of the given sample to
% be drawn.
% 
% PARAMETERS:
%   samples            - N samples drawn from the probability distribution
%   distr_fun          - function representing the probability distribution.
%                        The function should take two arguments (x, p) where x is 
%                        the argument to the distribution function and p the parameter(s)
%                        to be estimated
%   distr_fun_grad     - function representing the gradient of the distribution.
%                        It should take the same arguments as 'distr_fun'.
%   init_p             - initial guess for the parameter(s)
%   pmin               - minimum allowed value(s) for the parameter(s)
%   pmax               - maximum allowed value(s) for the parameter(s)
%
% RETURNS:
%   params  - the maximum-likelihood-estimation of the parameter(s), given
%             the samples

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

   mrstModule add optimization

   pmin = pmin(:);
   pmax = pmax(:);
   span = pmax - pmin;
   
   vertical = @(vec) vec(:);
   
   distr_fun_safe = @(x, p) max(eps, distr_fun(x, p));

   % NOTE: The reason the function works with the _logarithm_ of the probability
   % distribution rather than the probability distribution itself is better
   % numerical performance for distributions from the exponential family.  Such
   % distributions vanish so fast that their tails quickly become numerically equal
   % to zero, which may lead to problem when searching for optimal parameters.
   fun = @(u) vertical(sum(log(distr_fun_safe(samples, pmin + u .* span))));

   dfun = @(u) vertical(1./distr_fun_safe(samples, pmin + u(:) .* span))' * ...
               distr_fun_grad(vertical(samples), pmin + u(:) .* span);
   obj = @(u) deal(fun(u), (dfun(u)' .* span));

   start = (init_p(:) - pmin) ./ span;
   
   [v, u, history] = unitBoxBFGS(start, obj, 'plotEvolution', false); %#ok
   
   params = pmin + u .* span;
   
end
