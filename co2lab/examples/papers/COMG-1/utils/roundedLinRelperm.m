function fn = roundedLinRelperm(resSat, maxSat, epsilon)

   epsilon = min(epsilon, resSat/5);
   
   if epsilon == 0 
      % degenerate case (no residual saturation)
      fn = @(s) s;
   else 
      splinefn = compute_spline_function(resSat, maxSat, epsilon);
      fn = @(s) relperm(s, splinefn, resSat, maxSat, epsilon);
   end
end

% ----------------------------------------------------------------------------

function res = relperm(s, splinefn, resSat, maxSat, epsilon)

   curve_ix = (s >= (resSat - epsilon) & s <= (resSat + epsilon));
   lin_ix = (s > resSat + epsilon);
   
   res = 0*s;
   if sum(lin_ix) > 0
      res(lin_ix) = min((s(lin_ix) - resSat) .* maxSat ./ (maxSat - resSat), 1);
   end
   if sum(curve_ix) > 0
      res(curve_ix) = splinefn(s(curve_ix));
   end
   
   
end

% ----------------------------------------------------------------------------

function fun = compute_spline_function(resSat, maxSat, epsilon)
   
   p0 = 0;  % value of function at first point (resSat - epsilon)
   d0 = 0;  % derivative of function at first point (resSat - epsilon)
   p1 = epsilon * maxSat / (maxSat - resSat); % value of function at second point 
   d1 = maxSat / (1 - resSat);
   
   rhs = [p0; p1; d0; d1];
   
   xme = resSat - epsilon;
   xpe = resSat + epsilon;
   
   M = [xme^3,   xme^2, xme^1, 1; ...
        xpe^3,   xpe^2, xpe^1, 1; ...
        3*xme^2, 2*xme,     1, 0; ...
        3*xpe^2  2*xpe,     1, 0];
   
   coefs = M\rhs;
   a = coefs(1); b = coefs(2); c = coefs(3); d = coefs(4);

   fun = @(x) x.* ( x.* (a * x + b) + c) + d;
   
end

   