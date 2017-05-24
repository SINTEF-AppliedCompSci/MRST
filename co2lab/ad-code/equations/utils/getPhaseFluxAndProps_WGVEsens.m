function [vX, bX, mobX, rhoX, upc, dpX] = ...
   getPhaseFluxAndProps_WGVEsens(model, pW, pG, krX, T, gdz, phase, rs, temp,rhofac)%,permfac,porofac)
   % Function to compute phase fluxes, volume factors, mobilities,
   % densities, upstream indices and pressure gradients

   fluid = model.fluid;
   s     = model.operators;

   % Properties always computed using water pressure 
   switch upper(phase)
     case 'G'
       [bfun, mufun] = bind_temperature(fluid.bG, fluid.muG, temp);
       bX   = bfun(pW);
       rhoX = rhofac.*(bX .* fluid.rhoGS);
       mobX = krX ./ mufun(pW);
       pX   = pG;
     case 'W'
       [bfun, mufun] = bind_temperature(fluid.bW, fluid.muW, temp);
       bX   = bfun(pW);
       rhoX = bX .* (fluid.rhoWS + rs * fluid.rhoGS);
       mobX = krX ./ mufun(pW);
       pX   = pW;
     otherwise
       error('Indicated phase must be ''G'' (gas) or ''W'' (water)');
   end

   % rhoX on face, average of neighboring cells
   rhoXf = s.faceAvg(rhoX);

   % Compute pressure gradient, also taking into account the apparent gravity
   % component resulting from nonflat caprock geometry.
   dpX  = s.Grad(pX) - rhoXf .* gdz;

   % Identifying upstream side
   upc = double(dpX) <= 0;

   % Computing flux
   vX = -s.faceUpstr(upc, mobX) .* T .* dpX;
   if any (bX < 0)
      warning('Negative water compressibility present!');
   end
end

% ----------------------------------------------------------------------------

function [bfun, mufun] = bind_temperature(bfun, mufun, temp)

% A property function will be considered as dependent of temperature if it
% takes two arguments (not counting 'varargin')
   if (nargin(bfun) > 1) || (nargin(bfun) < -2)
      bfun = @(p) bfun(p, temp);
   end
   if (nargin(mufun) > 1) || (nargin(bfun) < -2)
      mufun = @(p) mufun(p, temp);
   end
end
