function Nc = computeCapillaryNumber(p, c, fluid, operators)

   s = operators;
   v = -s.T.*s.Grad(p);
   veloc = s.veloc;
   abs_veloc = 0;
   for i = 1 : numel(veloc)
      abs_veloc = abs_veloc + veloc{i}(v).^2;
   end
   abs_veloc = (abs_veloc).^(1/2);
   sigma = fluid.ift(c);
   Nc = abs_veloc./sigma;

end
