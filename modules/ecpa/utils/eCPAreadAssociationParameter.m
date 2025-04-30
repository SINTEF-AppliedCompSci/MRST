function ap = eCPAreadAssociationParameter(fluid, T)

na = fluid.names;
ncomp = numel(na);
ap = zeros(ncomp,ncomp);
for i = 1:ncomp
    for j = 1:ncomp
        if (strcmpi(na{i}, 'Water') && strcmpi(na{j}, 'CarbonDioxide')) ...
                || (strcmpi(na{j}, 'Water') && strcmpi(na{i}, 'CarbonDioxide'))
            ap(i,j) = 0.8061*(T/304.1282)-0.6904;
        elseif (strcmpi(na{i}, 'Water') && strcmpi(na{j}, 'HydrogenSulfide')) ...
                || (strcmpi(na{j}, 'Water') && strcmpi(na{i}, 'HydrogenSulfide'))
            ap(i,j) = 1.1472*(T/373.1)-0.6865;
        end
    end
end
end
