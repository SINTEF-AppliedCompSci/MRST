function state = updateExplicitAds(state, model, varargin)
% Explicit update of the concentration for the adsorption term
% We update the surfactant concentration ensuring mass conservation

   opt = struct('desorptionThreshold', -inf);
   opt = merge_options(opt, varargin{:});

   s = model.operators; % shortcut
   fluid = model.fluid; % shortcut

   [p, sW, c, cmax, ads, adsmax ] = model.getProps(state, 'pressure', 'water', 'surfactant', ...
                                                          'surfactantmax', 'ads', 'adsmax');

   new_ads = computeEffAds(c, opt.desorptionThreshold, fluid);

   % Multipliers for properties
   [pvMult, ~, ~] = getMultipliers(model.fluid, p, p);
   bW   = fluid.bW(p);
   poro = model.rock.poro;

   % Solve the equation
   %
   %   (s.pv/dt).*(pvMult.*bW.*sW.*(new_c - c) + fluid.rhoRSft.*((1-poro)./poro).*(new_ads - ads))= 0
   %
   % We remove values where sW is zero. At those points, the concentration is not defined and we
   % let it unchanged.
   %
   ind = sW > 1e-12;
   coef1 = fluid.rhoRSft.*((1-poro)./poro);
   coef2 = pvMult.*bW.*sW;
   coef = coef1./coef2;
   coef(~ind) = 0; % if sW=0, concentration is in fact not well-defined and we let it unchanged.
   new_c  = c - coef.*(new_ads - ads);

   negc = new_c < 0;
   % For the terms where the concentration becomes negative, we set the adsorption term so that
   % surfactant mass conservation is ensured, that is
   %
   %   (s.pv/dt).*(pvMult.*bW.*sW.*(0 - c) + fluid.rhoRSft.*((1-poro)./poro).*(new_ads - ads) )= 0
   %
   if any(negc)
      new_ads(negc) = ads(negc) + coef2(negc)./coef1(negc).*c(negc);
      new_c = max(0, new_c);
   end

   state = model.setProp(state, 'ads', new_ads);
   state = model.setProp(state, 'adsmax', max(new_ads, adsmax));
   state = model.setProp(state, 'surfactant', max(new_c, 0));
   state = model.setProp(state, 'surfactantmax', max(new_c, cmax));

end
