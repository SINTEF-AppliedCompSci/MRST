function dp = wellBoreFriction(v, rho, mu, D, L, roughness, flowtype, assumeTurbulent)
% Empricial model for well-bore friction

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

if nargin < 8
    assumeTurbulent = false;
end
if nargin < 7
    flowtype = 'massRate';
end
if numel(D) == 2
    [Di, Do] = deal(D(1), D(2));
else
    [Di, Do] = deal(0, D);
end

switch flowtype
    case 'velocity'
        % do nothing
    case 'volumeRate'
        v = v./(pi*((Do/2).^2 - (Di/2).^2));
    case 'massRate'
        v = v./(pi*rho.*((Do/2).^2 - (Di/2).^2));
    otherwise
        error(['Unknown flow type: ', flowtype]);
end

% small perturbation to prevent division by zero
isZero = value(v) == 0;
v(isZero) = eps; 

re = abs(rho.*v.*(Do-Di)./mu);

f  = (-3.6*log(6.9./re+(roughness./(3.7*(Do))).^(10/9))/log(10)).^(-2);

if ~assumeTurbulent % divide into laminar, intermediate and tubulent regions as in ECLIPSE
    [re1, re2] = deal(2000, 4000);
    lam = re <= re1;
    tur = re >= re2;
    int = ~or(lam, tur);
    f1 = 16/re1;
    f2 = (-3.6*log10(6.9./re2+(roughness./(3.7*(Do))).^(10/9))).^(-2);
    
    if any(lam)
        f(lam) = 16./re(lam);
    end
    
    if any(int) % interpolate between f1 f2
        if numel(f2) > nnz(int), f2 = f2(int); end
        f(int) = f1 + ((f2-f1)/(re2-re1)).*(re(int)-re1);
    end
end 

dp = -(2*sign(value(v)).*L./(Do-Di)).*(f.*rho.*v.^2);

% set dp-value for zero velocity segments, but keep derivative
if isa(dp, 'ADI')
    dp.val(isZero) = 0;
else
    dp(isZero) = 0;
end
end

