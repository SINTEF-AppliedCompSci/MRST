function state = updateExplicitAds(state0, state, model, dt)
% Explicit update of the concentration for the adsorption term
% It is important

   s = model.operators; % shortcut
   fluid = model.fluid; % shortcut


   [p, sW, c ] = model.getProps(state, 'pressure', 'water', 'surfactant', 'surfactantmax');

   [p0, cmax0, ads0, adsmax0] = model.getProps(state0, 'pressure', 'surfactantmax', 'ads', ...
                                                       'adsmax');

   ads = fluid.surfads(c);
   if fluid.adsInxSft == 2
      ads = max(ads, adsmax0);
   end


   % Multipliers for properties
   [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
   bW   = fluid.bW(p);
   poro = model.rock.poro;
   ads_term = fluid.rhoRSft.*((1-poro)./poro).*(ads - ads0);
   % Solve the equation
   %
   %   (s.pv/dt).*(pvMult.*bW.*sW.*(new_c - c) + ads_term )= 0
   %
   % We remove values where sW is zero. At those points, the concentration is not defined and we
   % let it unchanged.
   %
   ind = sW > 1e-12;
   coef = 1./(pvMult.*bW.*sW);
   coef(~ind) = 0;
   newc  = c - coef.*ads_term;

   negc = newc < 0;
   % For the terms where the concentration becomes negative, we set the adsorption term so that
   % surfactant mass conservation is ensured, that is
   %
   %   (s.pv/dt).*(pvMult.*bW.*sW.*(0 - c) + fluid.rhoRSft.*((1-poro)./poro).*(ads - ads0) )= 0
   %
   if any(negc)
      coef1 = fluid.rhoRSft.*((1-poro)./poro);
      coef2 = pvMult.*bW.*sW;
      ads(negc) = ads0(negc) + coef2(negc)./coef1(negc).*c(negc);
      c = max(0, newc);
   else
      c = newc;
   end

   state = model.setProp(state, 'ads', ads);
   adsmax = max(ads, adsmax0);
   state = model.setProp(state, 'adsmax', adsmax);
   state = model.setProp(state, 'surfactant', max(c, 0));
   state = model.setProp(state, 'surfactantmax', max(c, cmax0));


end
