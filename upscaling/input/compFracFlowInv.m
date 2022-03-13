function f_inv = compFracFlowInv(krw, kro, mu, swco, sowcr)
% only valid for S in [swco, sowcr]

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


ntab  = numel(krw);
f_inv = cell ([1, ntab]);



for i = 1:ntab
   s = linspace(swco(i), 1-sowcr(i), 100)';
   kr1 = krw{i}(s);
   kr2 = kro{i}(1-s);

   f_w = (kr1/mu(1))./(kr1/mu(1)+kr2/mu(2));

   [a,k,l] = unique(f_w, 'first');

   f_inv{i} = @(f) interp_f_inv(f_w(k), s(k), f, swco(i), 1-sowcr(i));
end





end

% Help function
%--------------------------------------------------------------------------
function sat = interp_f_inv(fw , s, f, s_min, s_max)

sat = interpTable(fw, s, f);
% not really necessary since fw is always in [0 1] and fw^-1(0) = swco and
% fw^-1(1) = sowcr:
% sat = min(max(sat, s_min), s_max);

end

