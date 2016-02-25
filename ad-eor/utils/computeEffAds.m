function ads = computeEffAds(c, adsmax, fluid)
% Compute effective adsorption, depending of desorption is present or not

% Treat differently case where no surfactant is present.  Physically, it does not make sense not to
% have fluid.surfads(0) = 0 and the input data should satisfied this identity However, we treat
% separatly the case c = 0 because, usually, adsorption increases very rapidly and this can create
% unnecessary convergence problem in the case where there is no surfactant.

   isSft = (double(c) > 0);
   ads = 0*c;
   ads(isSft) = fluid.surfads(c(isSft));
   if fluid.adsInxSft == 2
      ads = max(ads, adsmax);
   end

end
