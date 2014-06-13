function [yi, dyidxi, dyidvi] = interpRegPVT(T, xi, vi, flag, reginx)
compDer = (nargout>1);
nreg = numel(reginx);

yi = zeros(size(xi));
if compDer
    dyidxi = zeros(size(xi));
    dyidvi = zeros(size(xi));
end

tabSat = cellfun(@(x)x.data(x.pos(1:end-1),:), T, 'UniformOutput', false);

if nreg > 1
%     reginxSat  = reginx(flag,:);
%     reginxUSat = reginx(~flag,:);
    reginxSat  = cellfun(@(v) v(flag(v)), reginx, 'UniformOutput', false);
    reginxUSat = cellfun(@(v) v(~flag(v)), reginx, 'UniformOutput', false);
else
    reginxSat{1}  = ':';
    reginxUSat{1} = ':';
end

if ~compDer
    yi(flag)  = interpReg(tabSat, xi(flag), reginxSat);
    yi(~flag) = interp2DPVT(T, xi(~flag), vi(~flag), reginxUSat);
else
    [yi(flag), dyidxi(flag)] = interpReg(tabSat, xi(flag), reginxSat);
    [yi(~flag), dyidxi(~flag), dyidvi(~flag)] = interp2DPVT(T, xi(~flag), vi(~flag), reginxUSat);
end

end


