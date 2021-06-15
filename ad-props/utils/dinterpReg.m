function dyi = dinterpReg(T, xi, reginx)
% Interpolate table with multiple regions

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

nreg = numel(reginx);
if nreg > 1
    dyi = zeros(size(xi));
end
for k = 1:nreg
    if ischar(reginx{k}) && strcmp(reginx{k}, ':') %for improved eff seperate this case
        dyi = dinterpq1(T{k}(:,1), T{k}(:,2), xi);
    elseif ~isempty(reginx{k})
        dyi(reginx{k}) = dinterpq1(T{k}(:,1), T{k}(:,2), xi(reginx{k}));
    end
end
