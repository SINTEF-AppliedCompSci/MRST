function [krW, krO] = computeRelPermSft(sW, Nc, fluid)

   logNc = log(Nc)/log(10);
   % We cap logNc (as done in Eclipse)
   logNc = min(max(-20, logNc), 20);
   m = fluid.miscfact(logNc);

   sWcon    = fluid.sWcon;    % Residual water saturation   without surfactant
   sOres    = fluid.sOres;    % Residual oil saturation     without surfactant
   sWconSft = fluid.sWconSft; % Residual water saturation   with    surfactant
   sOresSft = fluid.sOresSft; % Residual oil saturation     with    surfactant

   % Interpolated water/oil residual saturations
   sNcWcon = m.*sWcon + (1 - m).*sWconSft;
   sNcOres = m.*sOres + (1 - m).*sOresSft;

   sNcEff = (sW - sNcWcon)./(1 - sNcWcon - sNcOres);

   % Rescaling of the saturation, without surfactant
   sNcWnoSft = (1 - sWcon - sOres).*sNcEff + sWcon;
   [krNcWnoSft, krNcOnoSft] = fluid.relPerm(sNcWnoSft);

   % Rescaling of the saturation, with surfactant
   sNcWSft =  (1 - sWconSft - sOresSft).*sNcEff + sWconSft;
   [krNcWSft, krNcOSft] = fluid.relPermSft(sNcWSft);

   [krW, krO] = fluid.relPerm(sW);

   krW = m.*krNcWnoSft + (1 - m).*krNcWSft;
   krO = m.*krNcOnoSft + (1 - m).*krNcOSft;

end
