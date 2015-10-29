function fn = roundedLinRelperm(resSat, epsilon)

   epsilon = min(epsilon, resSat/5);
   
   if epsilon == 0 
      % degenerate case
      fn = @(s) relperm(s, @(x) 0, 0, 0);
   else 
      splinefn = compute_spline_function(resSat, epsilon);
      fn = @(s) relperm(s, splinefn, resSat, epsilon);
   end
end

% ----------------------------------------------------------------------------

function res = relperm(s, splinefn, resSat, epsilon)

   curve_ix = (s >= (resSat - epsilon) & s <= (resSat + epsilon));
   lin_ix = (s > resSat + epsilon);
   
   res = 0*s;
   if sum(lin_ix) > 0
      res(lin_ix) = (s(lin_ix) - resSat) ./ (1 - resSat);
   end
   if sum(curve_ix) > 0
      res(curve_ix) = splinefn(s(curve_ix));
   end
   
   
end

% ----------------------------------------------------------------------------

function fun = compute_spline_function(resSat, epsilon)
   
   p0 = 0;  % value of function at first point (resSat - epsilon)
   d0 = 0;  % derivative of function at first point (resSat - epsilon)
   p1 = epsilon / (1-resSat); % value of function at second point (resSat + epsilon)
   d1 = 1 / (1-resSat);
   
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

   