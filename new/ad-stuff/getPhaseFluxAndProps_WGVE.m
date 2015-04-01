function [vX, bX, mobX, rhoX, upc, dpX] = ...
   getPhaseFluxAndProps_WGVE(model, pW, pG, krX, T, gdz, phase)
   
   fluid = model.fluid;
   s     = model.operators;
   
   % Properties always computed using water pressure @@ why?
   switch upper(phase)
     case 'G'
       bX   = fluid.bG(pW);
       rhoX = bX .* fluid.rhoGS;
       mobX = krX ./ fluid.muG(pW);
       pX   = pG;
     case 'W'
       bX   = fluid.bW(pW);
       rhoX = bX .* fluid.rhoWS;
       mobX = krX ./ fluid.muW(pW);
       pX   = pW;
     otherwise
       error('Indicated phase must be ''G'' (gas) or ''W'' (water)');
   end
   
   % rhoX on face, average of neighboring cells
   rhoXf = s.faceAvg(rhoX);
   
   % Compute pressure gradient, also taking into account the apparent gravity
   % component resulting from nonflat caprock geometry.
   dpX  = s.Grad(pX) - rhoXf .* gdz; % @@ minus or plus? (new operators!)
   
   % Identifying upstream side
   upc = double(dpX) <= 0;
   
   % Computing flux
   vX = -s.faceUpstr(upc, mobX) .* T .* dpX;
   if any (bX < 0)
      warning('Negative water compressibility present!');
   end
end