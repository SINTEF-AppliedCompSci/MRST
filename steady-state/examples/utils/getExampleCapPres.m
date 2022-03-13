function pcow = getExampleCapPres(poro, perm, swir, sor, varargin)
% Returns a capillary-pressure curves for each region. Returns a cell array
% of size nreg, where pcow{r} is a matrix of size opt.nsat x 2.

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
swir = opt.swir.*ones(nreg,1);
sor  = opt.sor.*ones(nreg,1);


% Compute capillary pressure

sW   = linspace(0, 1, nsat)';
cosa = cos(pi*opt.cangle/180); % cosine of contact angle
pcow = cell(1, nreg);

for r = 1:nreg
   pcow{r} = nan(nsat, 2);
   J = getJfunc(swir(r), sor(r), sW);
   sWsc = linspace(swir(r), 1-sor(r), nsat)';
   pcow{r}(:,1) = sWsc; % scaled to saturation limits
   pcow{r}(:,2) = J * surft * cosa * (poro(r)./perm(r))^(1/2); % pcow
end


end


function J = getJfunc(swir, sor, sW)
% Generates a Leverett J-function curve from an analytical expression

% Hardcoded parameters
jmult = 0.01;
delta = 0.05;

% J-function
sW    = sW(:);
J     = jmult*( 1./((sW-swir+delta).^2) - 1./((1-sor-sW+delta).^2) );

end


