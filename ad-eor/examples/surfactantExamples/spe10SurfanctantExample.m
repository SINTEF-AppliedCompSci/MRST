function spe10SurfanctantExample()
% Surfactant example for a SPE10 layer.

   try
      require ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui spe10
   catch
      mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui spe10
   end

   % Setup a TwoPhaseOilWaterModel with some standard fluid properties.
   [state, model, schedule] = setupSPE10_AD('layers')


   fluid = model.fluid;

   % Setup the fluid properties for our case

   % Setup the relative permeabilities using Corey model for the fluid without surfactant and
   % saturated with surfactant

   % Without surfactant
   n      = 3;   % Corey coefficient
   sWcon  = 0.2; % Residual water saturation
   sOres  = 0.2; % Residual oil saturation
   krWres = 0.6; % Endpoint relperm for water
   krOres = 0.5; % Endpoint relperm for oil

   krW = coreyPhaseRelpermAD(n, sWcon, krWres, sWcon + sOres);
   krO = coreyPhaseRelpermAD(n, sOres, krOres, sWcon + sOres);

   fluid.relPerm = @(sW) (relPerm(sW, krW, krO));
   fluid.sWcon = sWcon;
   fluid.sOres = sOres;

   % With surfactant
   n         = 1.5;
   sWconSft  = 0.05;
   sOresSft  = 0.05;
   krWresSft = 1;
   krOresSft = 1;

   krWSft = coreyPhaseRelpermAD(n, sWconSft, krWresSft, sWconSft + sOresSft);
   krOSft = coreyPhaseRelpermAD(n, sOresSft, krOresSft, sWconSft + sOresSft);

   fluid.relPermSft = @(sW) (relPerm(sW, krW, krO));
   fluid.sWconSft   = sWconSft;
   fluid.sOresSft   = sOresSft;

   % Reference pressure
   pRef = 234*barsa;

   % Compressibilities
   bW0       = 1./(1.012); % reference formation volume factor
   cW        = 4.28e-5/barsa; % compressibility coefficient
   fluid.bW  = @(p) bW0*exp((p - pRef)*cW);
   fluid.muW = @(p) 0*p + 0.48*centi*poise;
   bO0       = 1./(1.065); % reference formation volume factor
   cO        = 6.65e-5/barsa; % compressibility coefficient
   fluid.bO  = @(p) bO0*exp((p - pRef)*cO);
   fluid.muO = @(p) 0*p + 5*centi*poise;

   % Densities
   fluid.rhoWS = 962; % kg/m^3
   fluid.rhoOS = 1080;

   % Rock compressibility
   cR = 3e-5/barsa;
   fluid.cR = cR;
   fluid.pvMultR = @(p)(1 + cR.*(p-pRef));

   surfst = [ [   0;  30; 100], ...
              [0.61; 0.8;   1] ...
            ];
   surfst = extendTab(surfst); % extend to constant values.
   fluid.ift = @(c) interpReg(surfst, c);

   model = OilWaterPolymerModel(model.G, model.rock, fluid);

   injeIJ  = [56  10];
   prodIJ  = [ 5 211];

end

function [krW, krO] = relPerm(sW, krW, krO)
   krW = krW(sW);
   krO = krO(1 - sW);
end
