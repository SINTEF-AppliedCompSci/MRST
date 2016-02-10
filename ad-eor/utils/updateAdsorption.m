function state = updateAdsorption(state0, state, model)

   c = model.getProps(state, 'surfactant');
   adsmax0 = model.getProps(state, 'adsmax');

   ads = computeEffAds(c, adsmax0, model.fluid);
   adsmax = max(ads, adsmax0);
   state = model.setProp(state, 'ads', ads);
   state = model.setProp(state, 'adsmax', adsmax);

end


