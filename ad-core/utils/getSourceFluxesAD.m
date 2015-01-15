function [qRes, cells] = getSourceFluxesAD(model, mob, b, s, src)
    nPh = sum(model.getActivePhases);
    assert(size(src.sat, 2) == nPh);
    
    
    cells = src.cell;
    nsrc = numel(cells);
    
    inj = src.rate > 0;
    qRes = cell(nPh, 1);
    
    if any(~inj)
        totMob = mob{1};
        for i = 2:nPh
            totMob = totMob + mob{i};
        end
    end
    
    for i = 1:nPh
        q = double2ADI(zeros(nsrc, 1), mob{i});
        
        if any(inj)
            c = cells(inj);
            % Injection rates are given in reservoir conditions
            q(inj) = b{i}(c).*src.rate(inj).*src.sat(inj, i);
        end
        
        if any(~inj)
            c = cells(~inj);
            sc = s{i}(c);
            % Production rates are given in reservoir conditions. Use
            % mobilities ratios to ensure that immobile fluids are not
            % removed from the reservoir.
            f = mob{i}(c)./totMob(c);
            q(~inj) = b{i}(c).*f.*src.rate(~inj).*sc;
        end
        qRes{i} = q;
    end
end
