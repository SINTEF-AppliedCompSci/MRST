function [krW, krO] = computeRelPermSft(sW, Nc, fluid)

   logNc = log(Nc)/log(10);
   % We cap logNc (as done in Eclipse)
   logNc = min(max(-20, logNc), 20);
   m = fluid.miscfact(logNc);
   
   try
      set(0, 'currentFigure', 3);
   catch
      figure(3)
   end
   plot(m.val);
   title('m');
   axis([0, 100, 0, 1]);
   drawnow
   

   sWcon    = fluid.sWcon;    % Residual water saturation   without surfactant
   sOres    = fluid.sOres;    % Residual oil saturation     without surfactant
   sWconSft = fluid.sWconSft; % Residual water saturation   with    surfactant
   sOresSft = fluid.sOresSft; % Residual oil saturation     with    surfactant

   % Interpolated water/oil residual saturations
   sNcWcon = m.*sWconSft + (1 - m).*sWcon;
   sNcOres = m.*sOresSft + (1 - m).*sOres;

   sNcEff = (sW - sNcWcon)./(1 - sNcWcon - sNcOres);

   % Rescaling of the saturation, without surfactant
   sNcWnoSft = (1 - sWcon - sOres).*sNcEff + sWcon;
   [krNcWnoSft, krNcOnoSft] = fluid.relPerm(sNcWnoSft);

   % Rescaling of the saturation, with surfactant
   sNcWSft =  (1 - sWconSft - sOresSft).*sNcEff + sWconSft;
   [krNcWSft, krNcOSft] = fluid.relPermSft(sNcWSft);

   krW = m.*krNcWSft + (1 - m).*krNcWnoSft;
   krO = m.*krNcOSft + (1 - m).*krNcOnoSft;


end
