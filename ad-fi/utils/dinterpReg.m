function dyi = dinterpReg(T, xi, reginx)
nreg = numel(reginx);
if nreg > 1
    dyi = zeros(size(xi));
end
for k = 1:nreg
    if reginx{k} == ':' %for improved eff seperate this case
        dyi = dinterpTable(T{k}(:,1), T{k}(:,2), xi);
    elseif ~isempty(reginx{k})
        dyi(reginx{k}) = dinterpTable(T{k}(:,1), T{k}(:,2), xi(reginx{k}));
    end
end
