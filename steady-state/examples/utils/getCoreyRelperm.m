function getCoreyRelperm(swir, sor, n, krwmax)
%Undocumented Utility Function

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
   'nsat',    20  ... % Number of saturation points
   );
opt = merge_options(opt, varargin{:});

% Check input

nreg = numel(swir);
assert(numel(sor)==nreg, 'Wrong size of sor');

assert(numel(n)==nreg || numel(n)==1, 'Wrong size of n');
assert(numel(krwmax)==nreg || numel(krwmax)==1, 'Wrong size of krwmax');
n      = opt.n.*ones(nreg,1);
krwmax = opt.krwmax.*ones(nreg,1);



% Compute relative permeability

sW   = linspace(0, 1, opt.nsat)';
kr = cell(1, nreg);

for r = 1:nreg
   kr{r} = nan(nsat, 2);
   sWsc = swir(r) + sW./(1-swir(r)-sor(r));
   kr{r}(:,1) = sWsc; % scaled to saturation limits
   kr{r}(:,2) = J * surft * cosa * (poro(r)./perm(r))^(1/2); % pcow
end



end
