function fn = coreyPhaseRelpermAD(n, sr, kwm, sr_tot)
% Make function handle for AD style fluids
%
% n - exponent
% sr - residual saturation
% kwm - endpoint relperm
% sr_tot - sum of sr for all phases present

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
    
    if nargin < 1
        n = 1;
    end
    if nargin < 2
        sr = 0;
    end
    if nargin < 3
        kwm = 1;
    end
    if nargin < 4
        sr_tot = sr;
    end
%     den = 1 - sr_tot;
%     fn = @(s) kwm*max(min(((s - sr)./den), 1), 0).^n;
    fn = @(s, varargin) coreyRelperm(s, n, sr, kwm, sr_tot);
end

function kr = coreyRelperm(s, n, sr, kwm, sr_tot)
    den = 1 - sr_tot;
    sat = ((s - sr)./den);
    if isa(sat, 'ADI')
        sat.val = max(min(sat.val, 1), 0);
    else
        sat = max(min(sat, 1), 0);
    end
    
    kr = kwm*sat.^n;
end
