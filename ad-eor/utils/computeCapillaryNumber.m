function Nc = computeCapillaryNumber(p, c, fluid, operators)

   s = operators;
   v = -s.T.*s.Grad(p);
   veloc = s.veloc(v);
   abs_veloc = (sum(veloc.^2, 2)).^(1/2);
   sigma = fluid.ift(c);
   Nc = abs_veloc./sigma;

end
