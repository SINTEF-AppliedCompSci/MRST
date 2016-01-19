function [krW, krO] = computeRelPermSft(sW, Nc, fluid)

   logNc = log(Nc)/log(10);
   % We cap logNc (as done in Eclipse)
   logNc = min(max(-20, logNc), 20);
   m = fluid.miscfact(logNc);

   sWcon    = fluid.sWcon;
   sOres    = fluid.sOres;
   sWconSft = fluid.sWconSft;
   sOresSft = fluid.sOresSft;

   sNcWcon = m.*sWcon + (1 - m).*sWconSft;
   sNcOres = m.*sOres + (1 - m).*sOresSft;

   sNcWnosurf = ((1 - sNcWcon - sNcOres)./(1 - sWcon - sOres)).*(sW - sWcon) + sNcWcon;
   [krNcWnosurf, krNcOnosurf] = fluid.relPerm(sNcWnosurf);

   sNcWsurf = ((1 - sNcWcon - sNcOres)./(1 - sWconSft - sOresSft)).*(sW - sWconSft) + sNcWcon;
   [krNcWsurf, krNcOsurf] = fluid.relPermSft(sNcWsurf);

   krW = m.*krNcWnosurf + (1 - m).*krNcWsurf;
   krO = m.*krNcOnosurf + (1 - m).*krNcOsurf;

end
