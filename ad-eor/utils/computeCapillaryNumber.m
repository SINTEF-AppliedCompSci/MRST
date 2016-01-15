function [Nc, NcW] = computeCapillaryNumber(model, p, c)

   fluid = model.fluid;
   s = model.operators;
   v = -s.T.*s.Grad(p);
   poro =  s.pv./G.cells.volumes;
   poroFace = s.faceAvg(poro);
   faceA = G.faces.areas(s.internalConn);
   Kdp = v./(poroFace .* faceA);
   Kdp = abs(Kdp); % absolute value.
   
   sigma = fluid.ift(c);

   Nc = Kdp./sigma;

end
