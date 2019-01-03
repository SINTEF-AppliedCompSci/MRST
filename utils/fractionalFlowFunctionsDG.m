function f = fractionalFlowFunctionsDG(mob)

    nPh = numel(mob);
    mobT = @(s, sT, c) 0;
    for phNo = 1:nPh
        mobT = @(s, sT, c) mobT(s, sT, c) + mob{phNo}(s{phNo}, sT, c(:,phNo));
    end

    f = cell(nPh,1);
    for phNo = 1:nPh
        f{phNo} = @(s, sT, c) mob{phNo}(s{phNo}, sT, c(:,phNo))./mobT(s, sT, c);
    end
    
end