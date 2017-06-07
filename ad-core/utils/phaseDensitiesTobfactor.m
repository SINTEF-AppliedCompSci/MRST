function b = phaseDensitiesTobfactor(rho, rhoS, dissolved)
    numPh = numel(rhoS);
    b = rho;
    
    for i = 1:numPh
        factor = rhoS(i);
        if ~isempty(dissolved)
            for j = 1:numPh
                r_ph = dissolved{j}{i};
                if ~isempty(r_ph)
                    factor = factor + rhoS(j).*r_ph;
                end
            end
        end
        if iscell(b)
            b{i} = rho{i}./factor;
        else
            b(:, i) = rho(:, i)./factor;
        end
    end
end