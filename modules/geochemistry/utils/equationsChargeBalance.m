function [eqs, names, types] = equationsChargeBalance(model, state)
    
    chemsys = model.chemicalSystem;
    
    logSpecies  = model.getPropAsCell(state, 'logSpecies');
    species = cellfun(@(x) exp(x), logSpecies, 'UniformOutput', false);

    CV = chemsys.chargeVector;
    eInd = strcmpi('e-', chemsys.speciesNames);
    CV(1, eInd) = 0;

    CVp = max(CV, 0);
    CVn = max(-CV, 0);
    
    pos = 0;
    neg = 0;
    
    for k = 1 : chemsys.nC
        pos = pos + CVp(1, k).*species{k};
        neg = neg + CVn(1, k).*species{k};
    end

    eqs   = {log(pos) - log(neg)};
    names = {'charge balance equation'};
    types = {'cell'};
end

