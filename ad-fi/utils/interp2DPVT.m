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

        ixl = (bin==(tn-1));
        ixr = (bin==tn);
        if ~(ischar(reginx{k}) && reginx{k} == ':')
%             ixl = ixl.*reginx{k};
%             ixr = ixr.*reginx{k};
            ixl = reginx{k}(ixl);
            ixr = reginx{k}(ixr);        
        end
        
        
        yil(ixl) = interpTable(tab(:,1), tab(:,2), xi(ixl));
        yir(ixr) = interpTable(tab(:,1), tab(:,2), xi(ixr));
        if compDer
            dyidxil(ixl) = dinterpTable(tab(:,1), tab(:,2), xi(ixl));
            dyidxir(ixr) = dinterpTable(tab(:,1), tab(:,2), xi(ixr));
        end
    end
end
yi = yil.*w + yir.*(1-w);
if compDer
    dyidxi = dyidxil.*w + dyidxir.*(1-w);
    dyidvi = (yil-yir).*dwdvi;
end
end
