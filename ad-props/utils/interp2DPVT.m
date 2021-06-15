function [yi, dyidxi, dyidvi] = interp2DPVT(T, xi, vi, reginx)
compDer = (nargout>1);
nreg = numel(reginx);

w   = zeros(size(xi));
yil = zeros(size(xi));
yir = zeros(size(xi));
if compDer
    dyidxil = zeros(size(xi));
    dyidxir = zeros(size(xi));
    dwdvi   = zeros(size(xi));
end

for k = 1:nreg
    v = T{k}.key;  pos = T{k}.pos;
    lims = v; lims(1)=-inf; lims(end)=inf;
    vi_k = vi(reginx{k});
    [bin,bin] = histc(vi_k, lims);                                     %#ok
    w(reginx{k}) = (vi_k-v(bin))./(v(bin+1)-v(bin));
    if compDer
        dwdvi(reginx{k}) = 1./(v(bin+1)-v(bin));
    end
    for tn = 1:numel(v)
        lns = pos(tn):(pos(tn+1)-1);
        tab = T{k}.data(lns,:);
        
        if size(tab, 1) > 1 && tab(1, 1) > tab(2, 1)
            % If the data is decreasing rather than increasing, we reverse
            % the table.
            tab = tab(end:-1:1, :);
        end

        ixl = (bin==(tn-1));
        ixr = (bin==tn);
        if ~(ischar(reginx{k}) && strcmp(reginx{k}, ':'))
            ixl = reginx{k}(ixl);
            ixr = reginx{k}(ixr);        
        end
        
        
        yil(ixl) = fastInterpTable(tab(:,1), tab(:,2), xi(ixl));
        yir(ixr) = fastInterpTable(tab(:,1), tab(:,2), xi(ixr));
        if compDer
            dyidxil(ixl) = dinterpq1(tab(:,1), tab(:,2), xi(ixl));
            dyidxir(ixr) = dinterpq1(tab(:,1), tab(:,2), xi(ixr));
        end
    end
end
yi = yil.*w + yir.*(1-w);
if compDer
    dyidxi = dyidxil.*w + dyidxir.*(1-w);
    dyidvi = (yil-yir).*dwdvi;
end
end

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
