classdef UnivariateSplineFunction < handle
% Class representing a univariate (one parameter) spline function
% SYNOPSIS:
%   spline = UnivariateSplineFunction(order, knots, coefs)
%
% DESCRIPTION:
%
% REQUIRED PARAMETERS:
%
%   order  - Spline order (polynomial degree + 1)
%   knots -  Knot vector (vector), or number of internal intervals (integer).  If
%            number of intervals is specified, the knot vector will be
%            automatically generated as regularly spaced over the unit interval,
%            with k-multiple knots at endpoints.
%   coefs  - Spline control coefficients (control points).  If empty, all
%            control points will be set to zero by default.
%
% RETURNS:
%   Class instance.
%
% SEE ALSO:
%   BivariateSplineFunction
   
properties
   order % spline order (polynomial degree + 1)
   knots % knot vector (indicates parameter values where polynomial segments
         % are joined 
   coefs % coefficients for spline basis functions (B-spline functions)
end

methods (Access=public)

   % --------------------------------------------------------------------%
   function obj = UnivariateSplineFunction(order, knots, coefs)
   
      obj.order = order;
      if numel(knots) == 1
         obj.knots = UnivariateSplineFunction.makeRegularKnotvec(order, knots);
      else
         assert(numel(knots) >= 2 * order)
         obj.knots = knots;
      end
      
      if exist('coefs', 'var') && ~isempty(coefs)
         assert(numel(coefs) == numel(obj.knots) - order)
         obj.coefs = coefs(:)';
      else
         obj.coefs = zeros(1, numel(obj.knots) - order);
      end
   end

   % --------------------------------------------------------------------%      
   function vals = evaluate(self, param, deriv)
   % Evaluate spline or one of its derivatives.  
   %
   % INPUT: 
   % 
   % param  - vector of parameter values
   % 
   % deriv  - the derivative to evaluate (zero means evaluation of the spline
   %          function value)
   % 
   % RETURNS: 
   % 
   % vals - vector of spline (derivative) values for the specified parameters
      
      % computing matrix of all basis function values
      bfunvals = UnivariateSplineFunction.evaluateBfuns(self.order, ...
                                                        self.knots, ...
                                                        param(:), deriv);
      vals = bfunvals * self.coefs';
   end

   % --------------------------------------------------------------------%      
   function approximate(self, sample_params, sample_values, repar_spline, smoothing)
   % Adapt the spline function to approximate a given number of scattered input
   % data using penalized least squares.
   % 
   % INPUT:
   %
   % sample_params - parameter values for the input values to approximate
   %
   % sample_values - input values to approximate
   %
   % repar_spline  - if 'true', the parameter span of the spline (the knot
   %                 vector) will be shifted and re-scaled to fit the parameter
   %                 range of the input values/parameters.  Otherwise, it
   %                 will be set to span the unit interval, and input values
   %                 will be reparameterized to fit this interval before the
   %                 approximation algorithm is run.
   % 
   % smoothing     - The weight of the thin-plate spline smoothing term.  A
   %                 value of zero corresponds to no smoothing.  This may
   %                 lead to an indeterminate system if the number of degrees
   %                 of freedom (B-splines) is too high compared to the
   %                 number and distribution of the input data.
   % 
   % RETURNS:
   %
   % No return value.  This function operates on the spline object itself.
      
      range = [min(sample_params), max(sample_params)];
      if repar_spline
         % adapt spline parameterization to the range spanned by the sample
         % points
         self.knots = self.makeRegularKnotvec(self.order, self.numIntervals());
         self.knots = self.knots * diff(range) + range(1);
      else
         % reparametrizing input points to fit the unit interval
         self.knots = self.makeRegularKnotvec(self.order, self.numIntervals());
         sample_params = (sample_params - range(1)) / diff(range);
      end
      
      % setting up and solving normal equations
      F = UnivariateSplineFunction.evaluateBfuns(self.order, ...
                                                 self.knots, ...
                                                 sample_params, 0);
      rhs = F' * sample_values;
      FtF = F' * F;
      
      tot_span = self.knots(end) - self.knots(1);
      % scaling weight with span^3.  Explanation: The derivative scales inversely with
      % span.  The double derivative scales inversely with span squared.  The
      % product of two double derivatives thus scales inversely with span^4.
      % Integrating over the span thus yields something that scales inversely
      % with span^3,
      M = FtF + smoothing * (tot_span^3) * self.smoothnessMatrix(); 
      
      self.coefs = (M\rhs)'; 
   
   end
   
   % --------------------------------------------------------------------%      
   function res = numIntervals(self)
   % RETURNS: Nubmer of internal intervals in the knot vector
      res = numel(self.knots) - 2 * self.order + 1;
   end
   
   
   % --------------------------------------------------------------------%
   function omega = smoothnessMatrix(self)
   % RETURNS: The smoothness matrix corresponding to this spline.  This is the
   %          matrix that goes into the penalizing term of the approximate()
   %          function, and is constructed from a variational formulation
   %          intended to minimize the bending energy of the function,
   %          i.e. the integrated squared value of its second derivative.
         omega = self.innerProductMatrix(self.order, self.knots, 2); % second derivatives
   end
   
end

methods (Access=protected)

end

methods(Static)

   % --------------------------------------------------------------------%
   function kvec = makeRegularKnotvec(order, num_intervals)
      % constructing evenly-spaced knot vector over the unit interval with
      % k-multiple knots at end
      kvec = [zeros(1, order-1), ...
              linspace(0, 1, num_intervals+1), ...
              ones(1, order-1)];
   end
   
   % --------------------------------------------------------------------%
   function bvals = evaluateBfuns(order, knots, par, derivs)
      % evaluate all (derivatives of) B-spline basis functions for a set of
      % parameters
      if derivs > 0 
         tmp = UnivariateSplineFunction.evaluateBfuns(order-1, knots, par, derivs-1);
         w = knots(order:end) - knots(1:end-order+1);
         w(w==0) = 1; % these weights won't matter, but we modify them to
                      % avoid undefined values later
         w = (order-1)./w;
         tmp = w .* tmp;
         bvals = tmp(:, 1:end-1) - tmp(:, 2:end);
         return;
      end
               
      % Truncate parameters to lie in valid domain
      par(par<knots(1)) = knots(1);
      par(par>knots(end)) = knots(end);
         
      % determine correct interval for each parameter (start knot indices)
      
      s_ix = arrayfun(@(p) find(knots <= p, 1, 'last'), par, 'uniformoutput', false);
      s_ix = [s_ix{:}]';
      ix_max = numel(knots) - order;
      s_ix(s_ix > ix_max) = ix_max;
      
      % establish basis function values for order=1
      nint = numel(knots) - 1; % number of formal intervals (may be 0)
      bvals = accumarray([(1:numel(s_ix))', s_ix], 1, [numel(par), nint]);
      
      % iteratively compute basis values up to final order
      for k = 2:order
         % compute weights
         denoms = knots(k:end) - knots(1:end-k+1); % span [t_i, t_{i+k-1}]
         
         w1 = max(par - knots(1:end-k), 0);
         w1(w1>denoms(1:end-1)) = 0;
         
         w2 = max(knots(k+1:end) - par, 0); 
         w2(w2>denoms(2:end)) = 0;
         
         denoms(denoms==0) = 1;
         
         w1 = w1 ./ denoms(1:end-1);
         w2 = w2 ./ denoms(2:end);

         bvals = bsxfun(@times, w1, bvals(:,1:end-1)) + ...
                 bsxfun(@times, w2, bvals(:,2:end));
      end
      
   end
   
   % ----------------------------------------------------------------------------
   function gq = gaussQuad3(order, knots, deriv)
      % For each interval, evaluate each of its nonzero B-spline basis
      % functions at the gauss-Legendre quadrature points (3-point version)
      
      N = numel(knots) - order; % number of B-spline basis functions
      gq = zeros(3*(order), N);
      
      for i = 1:N
         
         % compute quadrature points for basis spline i, extending from knot
         % t_i to t_{i+order}.
         im = 0.5 * (knots(i+1:i+order) + knots(i:i+order-1));  % midpoints
         iw = 0.5 * (knots(i+1:i+order) - knots(i:i+order-1));  % semi-widths
         
         par = [im - iw*sqrt(3/5); im; im + iw*sqrt(3/5)];
         par = par(:);
         
         % @@ Ineffective.  We evaluate for all basis function where we
         % only need the result for one.
         tmp = UnivariateSplineFunction.evaluateBfuns(order, knots, par, deriv);
         
         gq(:,i) = tmp(:, i);
      end
      
   end
   
   % ----------------------------------------------------------------------------
   function omega = innerProductMatrix(order, knots, deriv)
      % Returns the matrix with inner products of (B-spline) basis functions,
      % or of a given derivative thereof.
      
      assert(order - deriv > 0) % otherwise, the inner product will be zero
      assert(order - deriv <= 5) % otherwise, we are using too few quad. points
      
      % aligning computing quadrature points and aligning them row-wise in a
      % global matrix
      qpts = UnivariateSplineFunction.gaussQuad3(order, knots, deriv);
      nb = size(qpts, 2); % number of basis functions
      nq = size(qpts, 1); % number of quadrature points
      
      ix = repmat(1:nb, nq, 1);      
      cstart = 1:3:3*nb;
      jx = mcolon(cstart, cstart + nq - 1); 
      
      % one row per basis function, three columns per interval, representing
      % the three quadrature point evaluations for that basis function on
      % that interval
      qmat = accumarray([ix(:), jx(:)], qpts(:)');
      
      % truncating away the columns of omega that correspond to intervals
      % outside the main parameter range
      qmat = qmat(:, 3*(order-1)+1:end-3*(order-1));

      % spans of each interval 
      spans = knots(2:end) - knots(1:end-1);
      spans = spans(order:end-order+1);
      
      % weights
      w   = [5/9; 8/9; 5/9]; % quadrature weights
      sw = w * spans; % span multiplied by weight
      
      % computing diagonals
      d = zeros(nb, order);

      % main diagonal (square integral of basis function)
      % we divide the result by two, since we will add the resulting matrix
      % to its transpose later
      tmp = bsxfun(@times, qmat.^2, sw(:)');
      d(:,1) = 0.5 * (0.5 * sum(tmp, 2));
      
      % side diagonals (products between neighboring basis functions)
      for d_ix = 1:order-1
         tmp = qmat(1:end-d_ix, :) .* qmat(d_ix+1:end, :);
         tmp = bsxfun(@times, tmp, sw(:)');
         d(d_ix+1:end, d_ix+1) = 0.5 * sum(tmp, 2);
      end
      omega = spdiags(d, 0:order-1, nb, nb);
      omega = omega + omega';
   end
   
end

end

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
