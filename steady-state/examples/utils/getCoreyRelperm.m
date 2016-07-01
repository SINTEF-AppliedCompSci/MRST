function getCoreyRelperm(swir, sor, n, krwmax)
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
