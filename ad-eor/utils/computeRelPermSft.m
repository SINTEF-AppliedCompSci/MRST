function [krW, krO] = computeRelPermSft(fluid, sW, Nc)

   m = fluid.surfcapd(Nc);

   sWcon     = fluid.sWcon;
   sOres     = fluid.sOres;
   sWconSurf = fluid.sWconSurf;
   sOresSurf = fluid.sOresSurf;

   sNcWcon = m.*sWcon + (1 - m).*sWconSurf;
   sNcOres = m.*sOres + (1 - m).*sOresSurf;

   sNcWnosurf = ((1 - sNcWcon - sNcOres)./(1 - sWcon - sOres)).*(sW - sWcon) + sNcWcon;
   krNcWnosurf = fluid.krW(sNcWnosurf);
   krNcOnosurf = fluid.krO(sNcWnosurf);

   sNcWsurf = ((1 - sNcWcon - sNcOres)./(1 - sWconSurf - sOresSurf)).*(sW - sWconSurf) + sNcWcon;
   krNcWsurf = fluid.krW(sNcWsurf);
   krNcOsurf = fluid.krO(sNcWsurf);

   krW = m.*krNcWnosurf + (1 - m).*krNcWsurf;
   krO = m.*krNcOnosurf + (1 - m).*krNcOsurf;

end
