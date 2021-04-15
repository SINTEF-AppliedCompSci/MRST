function WP = modelAverageWellPairDiagnostics(d, modselIx)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

assert(numel(modselIx)>1);
if nargin < 2
    modselIx = (1:numel(d.Data))';
end

nM = numel(modselIx);

w  = 1/nM;

for i = 2:(numel(modselIx))
    diagn{i-1} = d.Data{modselIx(i)}.diagnostics;
end
WP    = diagn{1}.WP;
for k = 1:(numel(modselIx)-1)
    z = k>1;
    for ni = 1:numel(WP.inj)
        WP.inj(ni).alloc  = z*WP.inj(ni).alloc  + w*diagn{k}.WP.inj(ni).alloc;
        WP.inj(ni).ralloc = z*WP.inj(ni).ralloc + w*diagn{k}.WP.inj(ni).ralloc;
    end
    
    for np = 1:numel(WP.prod)
        WP.prod(np).alloc  = z*WP.prod(np).alloc  + w*diagn{k}.WP.prod(np).alloc;
        WP.prod(np).ralloc = z*WP.prod(np).ralloc + w*diagn{k}.WP.prod(np).ralloc;
    end
    
    WP.vols = z*WP.vols + w*diagn{k}.WP.vols;
end
