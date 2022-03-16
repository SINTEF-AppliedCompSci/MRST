function pdf = estimatePDF(mean, median, std, scale)
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

    gamma3 = 3*(mean-median)/std;
    gamma3 = sign(gamma3)*min(0.99, abs(gamma3));
%     gamma3 = -0.8
    a     = abs(gamma3)^(2/3);
    delta = sign(gamma3).*sqrt(pi/2*a/(a+((4-pi)/2)^(2/3)));
    alpha = delta/sqrt(1-delta^2);
    
    std = std/sqrt((1 - 2*delta^2/pi));
    mean = mean - std*delta*sqrt(2/pi);
    
    xi = @(x) (x-mean)./std;
    
    pdf = @(x) 2/(std*sqrt(2*pi))*exp(-0.5*xi(x).^2).*0.5.*(1 + erf(alpha*xi(x)/sqrt(2))).*scale;
    
end
