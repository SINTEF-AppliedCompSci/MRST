function ads = computeEffAds(c, adsmax, fluid)
   % Compute effective adsorption, depending of desorption is present or not
   ads = fluid.surfads(c);
   if fluid.adsInxSft == 2
      ads = max(ads, adsmax);
   end
end
