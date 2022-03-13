classdef BivariateSplineFunction < handle
% Class representing a bivariate (two parameter) tensor product spline function
%
% SYNOPSIS:
%   obj = BivariateSpline(orders, knots_u, knots_v, coefs)
%
% DESCRIPTION:
%
% PARAMETERS:
%   orders  - [order_u, order_v] - spline orders in each parameter
%   knots_u - Knot vector (vector) or number of internal intervals (integer) for
%             the u-parameter.  If number of intervals is specified, the knot
%             vector will be automatically generated as regularly spaced over
%             the unit interval, with k-multiple knots at the endpoints.
%   knots_v - Same description as above, but for the v-parameter
%   coefs   - Spline control coefficients (control points).  If empty, all
%             control points will be set to zero by default. 
%
% RETURNS:
%   Class instance.
%
%
% SEE ALSO:
%   UnivariateSplineFunction

   
properties
   orders   % [order_u, order_v]
   knots_u  % knot vector in u-parameter
   knots_v  % knot vector in v-parameter
   coefs    % coefficient grid
end

methods (Access=public)
   % --------------------------------------------------------------------%
   function obj = BivariateSplineFunction(orders, knots_u, knots_v, coefs)
      assert(numel(orders) == 2)
      obj.orders = orders;
      obj.knots_u = BivariateSplineFunction.init_kvec(orders(1), knots_u);
      obj.knots_v = BivariateSplineFunction.init_kvec(orders(2), knots_v);
      if exist('coefs', 'var') && ~isempty(coefs)
         assert(size(coefs, 1) == numel(obj.knots_u) - orders(1));
         assert(size(coefs, 2) == numel(obj.knots_v) - orders(2));
         obj.coefs = coefs(:);
      else
         obj.coefs = zeros(numel(obj.knots_u) - orders(1), ...
                           numel(obj.knots_v) - orders(2));
      end
      
   end

   % --------------------------------------------------------------------%
   function d = domain(self)
   % Return parameter domain on the form [umin, umax, vmin, vmax]
      d = [self.knots_u(self.orders(1)), ...
           self.knots_u(end - self.orders(1) + 1), ...
           self.knots_v(self.orders(2)), ...
           self.knots_v(end - self.orders(2) + 1)];
   end
                  
   % --------------------------------------------------------------------%
   function vals = evaluate(self, params, deriv1, deriv2)
   % Evaluate spline or one of its derivatives.  
   %
   % INPUT: 
   % 
   % params  - 2-column matrix with parameter values.  First column
   %           represents u-parameter values, second column v-parameter
   %           values. 
   % 
   % deriv1  - the derivative to evaluate (zero means evaluation of the spline
   %           function value) for the u-parameter
   % deriv2  - the derivative to evaluate (zero means evaluation of the spline
   %           function value) for the v-parameter
   % 
   % RETURNS: 
   % 
   % vals - vector of spline (derivative) values for the specified parameters
      
      if ~exist('deriv1', 'var')
         deriv1 = 0;
      end
      if ~exist('deriv2', 'var')
         deriv2 = 0;
      end
      
      bu = UnivariateSplineFunction.evaluateBfuns(self.orders(1), ...
                                                  self.knots_u, ...
                                                  params(:,1), deriv1);
      bv = UnivariateSplineFunction.evaluateBfuns(self.orders(2), ...
                                                  self.knots_v, ...
                                                  params(:,2), deriv2);
      
      % @@ check if the following loop can be optimized
      vals = zeros(size(params, 1), 1);
      for i = 1:size(params, 1)
         vals(i) = kron(bv(i,:), bu(i,:)) * self.coefs(:);
      end
      
   end
   
   % --------------------------------------------------------------------%
   function approximate(self, sample_params, sample_values, repar_spline, smoothing)
   % Adapt the spline function to approximate a given number of scattered input
   % data using penalized least squares.
   % 
   % INPUT:
   %
   % sample_params - parameter values for the input values to approximate
   %                 (two-column vector, one column for u-parameter, one for
   %                 v-parameter) 
   %
   % sample_values - input values to approximate
   %
   % repar_spline - 'none', 'data', 'unit'.  'none' means no
   %                 reparametrization (user responsible to ensure that
   %                 sample values all fall within the valid range).  'data'
   %                 means that the knot vectors will be shifted and rescaled
   %                 to fit the parameter range of the input data.  'unit'
   %                 means that the knot vectors will be set to span the unit
   %                 square (k-regular), and the input values will be
   %                 reparametrized to also fit in the unit square.
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
      
      assert(size(sample_params, 1) == size(sample_values,1));
      range_u = [min(sample_params(:,1)), max(sample_params(:,1))];
      range_v = [min(sample_params(:,2)), max(sample_params(:,2))];
      
      if repar_spline == 'data'
         % adapt spline parameterization to the range spanned by the sample
         % points
         self.knots_u = self.kvec_to_range(self.orders(1), ...
                                           self.numIntervals(1), range_u);
         self.knots_v = self.kvec_to_range(self.orders(2), ...
                                           self.numIntervals(2), range_v);
      elseif repar_spline == 'unit'
         self.knots_u = self.init_kvec(self.orders(1), self.numIntervals(1));
         self.knots_v = self.init_kvec(self.orders(2), self.numIntervals(2));
         sample_params(:,1) = (sample_params(:,1) - range_u(1)) / diff(range_u);
         sample_params(:,2) = (sample_params(:,2) - range_v(1)) / diff(range_v);
      elseif repar_spline ~= 'none'
         error('invalid option')
      end
      
      % setting up and solving normal equations
      n = numel(sample_values); % number of samples
      Fu = UnivariateSplineFunction.evaluateBfuns(self.orders(1), ...
                                                  self.knots_u, ...
                                                  sample_params(:,1), 0);
      Fv = UnivariateSplineFunction.evaluateBfuns(self.orders(2), ...
                                                  self.knots_v, ...
                                                  sample_params(:,2), 0);
      
      ix_i = [];
      ix_j = [];
      vals = [];
      for r = 1:n
         tmp = kron(Fv(r,:), Fu(r,:));
         ix = find(tmp);
         ix_i = [ix_i; r * ones(numel(ix), 1)];
         ix_j = [ix_j; ix(:)];
         vals = [vals; (tmp(ix))'];
      end
      F = accumarray([ix_i, ix_j], vals, [n, numel(self.coefs)], [], 0, true);
      
      % Double derivatives scale inversely with parameter span squared.
      % Integrating over span in two directions thus yields something that
      % scale inversely with (span)^3
      span_u = self.knots_u(end) - self.knots_u(1);
      span_v = self.knots_v(end) - self.knots_v(1);
      fac = max(span_u, span_v)^3; 
      
      rhs = F' * sample_values;
      FtF = F' * F;
      M = FtF + smoothing * fac * self.smoothnessMatrix();
      
      self.coefs = reshape((M\rhs), size(self.coefs, 1), size(self.coefs, 2));
   end
   
   % --------------------------------------------------------------------%
   function res = numIntervals(self, ix)
      % RETURNS: Nubmer of internal intervals in the u (ix=1) or v (ix=2)
      % knot vector.
      assert(ix == 1 || ix == 2); % either u or v direction

      if ix == 1
         res = numel(self.knots_u) - 2 * self.orders(1) + 1;
      else
         res = numel(self.knots_v) - 2 * self.orders(2) + 1;
      end
   end
   
   % --------------------------------------------------------------------%
   function omega = smoothnessMatrix(self)
   % RETURNS: The smoothness matrix corresponding to this spline function.  This is
   % the matrix that goes into the penalizing term of the approximate()
   % function, and is constructed from a variational formulation intended to
   % minimize the bending energy of the function, i.e. the integrated squared
   % value of its second derivative.

      % shortcuts
      [ou, ov] = deal(self.orders(1), self.orders(2));
      [ku, kv] = deal(self.knots_u, self.knots_v);
      
      % matrix from differentiating the (d2/dx2 f)^2 part
      m1 = kron(UnivariateSplineFunction.innerProductMatrix(ov, kv, 0), ...
                UnivariateSplineFunction.innerProductMatrix(ou, ku, 2));
      
      % matrix from differentiating the (d/dx d/dy f)^2 part
      m2 = kron(UnivariateSplineFunction.innerProductMatrix(ov, kv, 1), ...
                UnivariateSplineFunction.innerProductMatrix(ou, ku, 1));
      
      % matrix from differentiating the (d2/dy2 f)^2 part
      m3 = kron(UnivariateSplineFunction.innerProductMatrix(ov, kv, 2), ...
                UnivariateSplineFunction.innerProductMatrix(ou, ku, 0));
      
      omega = m1 + 2 * m2 + m3;
   end      
end

methods (Access=protected)
   
end

methods (Static)
   % --------------------------------------------------------------------%
   function kvec = kvec_to_range(order, num_intevals, range)
      kvec = UnivariateSplineFunction.makeRegularKnotvec(order, num_intevals);
      kvec = kvec * diff(range) + range(1);
   end

   % --------------------------------------------------------------------%
   function kvec = init_kvec(order, input)
      if numel(input) == 1
         kvec = UnivariateSplineFunction.makeRegularKnotvec(order, input);
      else
         assert(numel(input) >= 2 * order)
         kvec = input(:)';
      end
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
