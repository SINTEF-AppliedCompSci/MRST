function p = getExampleFluidProps(varargin)
% Returns a structure with fluid properties for use with the examples.
% 
% A fluid can be created from the returned structe by calling the function 
% initADIFluidOWPolymer.
% 
opt = struct(...
	'satnum', [],  ...
    'swir',   0.1, ... % Irreducible water saturation
    'sor',    0.1, ... % Residual oil saturation
    'krexp',  2,   ... % Relperm exponent
    'krWmax', 0.5  ... % Maximum value of water relperm
    );
opt = merge_options(opt, varargin{:});

% Region map
p.satnum = satnum;

% Oil and Water Properties ------------------------------------------------

% Residual Water Saturation
p.swir = 0.1;
p.sor  = 0.1;

% Viscosity (cP)
% Function of pressure (psia)
p.muW = 1.0*centi*poise;
p.muO = 1.5*centi*poise;

% Density (kg/m^3)
% Function of pressure (psia)
p.rhoW = 1014*kilogram/meter^3;
p.rhoO = 850*kilogram/meter^3;

% Relperm (unitless)
% Function of phase saturation (unitless)
sW    = linspace(0, 1, 20)';
sWlim = linspace(p.swir, 1-p.sor, 20)';
p.krW = [   sWlim   opt.krWmax.*(sW.^opt.krn)];
sW    = flipud(sW);
sWlim = flipud(sWlim);
p.krO = [(1-sWlim) (1-sW).^opt.krn];

% Capillary Pressure
% Parameters for j-function convertion
poro = 0.1;
perm = 100.*milli*darcy;
p.pcOW = pcowCurve(poro, perm, 'swir', p.swir, 'sor', p.sor, ...
   'nsat', numel(sW));

% Formation Volume Factor (unitless)
% Function of pressure (psia)
p.BW = 1;
p.BO = 1;
% p.BW = {'fvf', [200*barsa 1.0 4.28e-5*(1/barsa)]}; % pref, fvfref, comp
% p.BO = {'fvf', [200*barsa 1.0 6.65e-5*(1/barsa)]}; % pref, fvfref, comp

% Pore Volume Multiplier (unitless)
p.pvMultR = [];



% Polymer Properties ------------------------------------------------------

% Adsorption (kg/kg)
% Function of polymer concentration (kg/m^3)
p.ads = [ [ 0; 2.3;  3;  4].*kilogram/meter^3, ...
          [ 0;  20; 20; 20].*(milli*gram)/(kilo*gram) ];

% Maximum Adsorption (kg/kg)
p.adsMax = 20.*(milli*gram)/(kilo*gram);

% Set Desorption
% adsInx=1: desorption
% adsInx=2: no desorption
p.adsInx = 1;

% Viscosity Multiplier (unitless)
% Function of polymer concentration (kg/m^3)
p.muWMult = [ [ 0;  3;  4].*kilogram/meter^3, ...
              [ 1; 40; 40] ];

% Dead Pore Space / Inaccesible Pore Volume (unitless)
p.dps = 0.0;

% Residual Resistance Factor, RRF (unitless)
p.rrf = 1.3;

% Rock Density (kg/m^3)
p.rhoR = 2000.*kilogram/meter^3;

% Todd-Longstaff Mixing Parameter, Omega (unitless)
p.mixPar = 1;

% Maximum Polymer Consentration (kg/m^3)
p.cmax = 4;

end



%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------




