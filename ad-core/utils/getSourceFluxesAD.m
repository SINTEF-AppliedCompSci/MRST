function [qRes, cells] = getSourceFluxesAD(model, mob, b, s, src)
    nPh = sum(model.getActivePhases);
    assert(size(src.sat, 2) == nPh);
    
    
    cells = src.cell;
    nsrc = numel(cells);
    
    inj = src.rate > 0;
    qRes = cell(nPh, 1);
    for i = 1:nPh
        q = double2ADI(zeros(nsrc, 1), mob{i});
        
        if any(inj)
            q(inj) = src.rate(inj).*src.sat(inj, i);
        end
        
        if any(~inj)
            c = cells(~inj);
            sc = s{i}(c);
            q(~inj) = b{i}(c).*src.rate(~inj).*sc;
        end
        qRes{i} = q;
    end
end