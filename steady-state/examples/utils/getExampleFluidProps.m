function p = getExampleFluidProps(rock, varargin)
% Returns a structure with fluid properties for use with the examples.
% 
% A fluid can be created from the returned structe by calling the function 
% initADIFluidOWPolymer.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

opt = struct(...
	'satnum',  [],   ...
    'nsat',    20,   ... % Number of saturaton points
    'swir',    0.1,  ... % Irreducible water saturation
    'sor',     0.1,  ... % Residual oil saturation
    'krn',     2,    ... % Relperm exponent    
    'krWmax',  0.5,  ... % Maximum value of water relperm
    'polymer', false ...
    );
opt = merge_options(opt, varargin{:});

% Region map
p.satnum = opt.satnum;
nreg = max(p.satnum);
assert(numel(unique(p.satnum)) == nreg, ...
    'Satnum regions must be numbered from 1 through nreg');

% Ensure values have the correct length. If only a single value is given,
% the same value is assumed for all regions.
opt.swir   = opt.swir(:).*ones(nreg,1);
opt.sor    = opt.sor(:).*ones(nreg,1);
opt.krn    = opt.krn(:).*ones(nreg,1);
opt.krWmax = opt.krWmax(:).*ones(nreg,1);


% Oil and Water Properties ------------------------------------------------

% Residual Water Saturation
p.swir = opt.swir;
p.sor  = opt.sor;

% Viscosity (cP)
% Function of pressure (psia)
p.muW = 1.0*centi*poise;
p.muO = 1.5*centi*poise;

% Density (kg/m^3)
% Function of pressure (psia)
p.rhoW = 1014*kilogram/meter^3;
p.rhoO = 850*kilogram/meter^3;

% Corey-type Relperms (unitless)
% Function of phase saturation (unitless)
sW  = linspace(0, 1, opt.nsat)';
krW = cell(1,nreg);
krO = cell(1,nreg);
for r = 1:nreg
    sWlim  = linspace(p.swir(r), 1-p.sor(r), opt.nsat)';
    krW{r} = [   sWlim   opt.krWmax(r).*(sW.^opt.krn(r))];
    krO{r} = flipud([(1-sWlim) (1-sW).^opt.krn(r)]);
end
if r==1
    krW = krW{1};
    krO = krO{1};
end
p.krW = krW;
p.krO = krO;

% Capillary Pressure
% Parameters for j-function convertion
poro = nan(nreg,1);
perm = nan(nreg,1);
for r = 1:nreg
    poro(r) = mean(rock.poro(p.satnum==r));
    perm(r) = 10.^mean(log10(rock.perm(p.satnum==r,1)));
end
p.pcOW = getCapPres(poro, perm, p.swir, p.sor, 'nsat', opt.nsat);

% Formation Volume Factor (unitless)
% Function of pressure (psia)
p.BW = 1;
p.BO = 1;
% p.BW = {'fvf', [200*barsa 1.0 4.28e-5*(1/barsa)]}; % pref, fvfref, comp
% p.BO = {'fvf', [200*barsa 1.0 6.65e-5*(1/barsa)]}; % pref, fvfref, comp

% Pore Volume Multiplier (unitless)
p.pvMultR = [];



% Polymer Properties ------------------------------------------------------
%
% There are polymer properties for one or two regions.
% 

if opt.polymer

    % Adsorption (kg/kg)
    % Function of polymer concentration (kg/m^3)
    p.ads = {[ [ 0; 2.3;  3;  4].*kilogram/meter^3, ...
               [ 0;  20; 20; 20].*(milli*gram)/(kilo*gram) ], ...
             [ [ 0;   2;  3;  4].*kilogram/meter^3, ...
               [ 0;  70; 70; 70].*(milli*gram)/(kilo*gram) ] };
    
    % Maximum Adsorption (kg/kg)
    p.adsMax = [40; 140].*(milli*gram)/(kilo*gram);
    
    % Set Desorption
    % adsInx=1: desorption
    % adsInx=2: no desorption
    p.adsInx = 1;

    % Viscosity Multiplier (unitless)
    % Function of polymer concentration (kg/m^3)
    p.muWMult = {[ [ 0;   3;  4].*kilogram/meter^3, ...
                   [ 1;  40; 40] ], ...
                 [ [ 0;   3;  4].*kilogram/meter^3, ...
                   [ 1;  30; 30] ] };

    % Dead Pore Space / Inaccesible Pore Volume (unitless)
    p.dps = 0.0;

    % Residual Resistance Factor, RRF (unitless)
    p.rrf = [1.1; 1.3];

    % Rock Density (kg/m^3)
    p.rhoR = 2000.*kilogram/meter^3;

    % Todd-Longstaff Mixing Parameter, Omega (unitless)
    p.mixPar = 1;

    % Maximum Polymer Consentration (kg/m^3)
    p.cmax = 4;
    
    % If there is only a single region, we remove the properties of the
    % second region
    if (nreg==1)
        p.ads     = p.ads{1};
        p.adsMax  = p.adsMax(1);
        p.muWMult = p.muWMult{1};
    end
    
end

end



%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------


function pcow = getCapPres(poro, perm, swir, sor, varargin)
% Returns a capillary-pressure curves for each region. Returns a cell array
% of size nreg, where pcow{r} is a matrix of size opt.nsat x 2.
opt = struct(...
   'surft',   22*dyne/(centi*meter), ... % Surface tension 
   'cangle',  45, ... % Contact angle (degrees)
   'nsat',    20  ... % Number of saturation points
   );
opt = merge_options(opt, varargin{:});


% Check input

surft = opt.surft;
assert(isscalar(surft), 'Surface tension must be scalar');

nsat = opt.nsat;

nreg = numel(perm);
assert(numel(poro)==nreg, ...
   'Permeability and porosity must be of same size');

assert(numel(swir)==nreg || numel(swir)==1, 'Wrong size of swir');
assert(numel(sor)==nreg || numel(sor)==1, 'Wrong size of sor');
swir = swir.*ones(nreg,1);
sor  = sor.*ones(nreg,1);


% Compute capillary pressure

cosa = cos(pi*opt.cangle/180); % cosine of contact angle
pcow = cell(1, nreg);
J    = getJfunc(nsat);

for r = 1:nreg
   pcow{r} = nan(nsat, 2);
   sWsc = linspace(swir(r), 1-sor(r), nsat)';
   pcow{r}(:,1) = sWsc; % scaled to saturation limits
   pcow{r}(:,2) = J * surft * cosa * (poro(r)./perm(r))^(1/2); % pcow
end


end


function J = getJfunc(npoints)
% Generates a Leverett J-function curve from an analytical expression

% Hardcoded parameters
delta = 0.05; % Determines the "size of the asymptotes"
Jmax  = 0.9;
Jmin  = -1.8;

% J-function
sW = linspace(delta,1-delta,npoints)';
JL = 1./((sW).^2);   JL = Jmax.*( JL./max(JL) );
JR = 1./((1-sW).^2); JR = Jmin.*( JR./max(JR) );
J  = JL + JR;

end




