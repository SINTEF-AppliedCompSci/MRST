function [eqs, names, types] = equationsChargeBalance(model, state)
    
    chemsys = model.chemicalSystem;
    
    logSpecies  = model.getPropAsCell(state, 'logSpecies');
    logElements = model.getPropAsCell(state, 'logElements');
    if chemsys.nLC > 0
        combinationComponents = model.getPropAsCell(state, 'combinationComponents');
    end
    if chemsys.nG > 0
        logPartialPressures = model.getPropAsCell(state, 'logPartialPressures');
    end
    if chemsys.nS > 0    
        logSaturationIndicies = model.getPropAsCell(state, 'logSaturationIndicies');
    end
    if chemsys.nP > 0    
        logSurfaceActivityCoefficients = model.getPropAsCell(state, 'logSurfaceActivityCoefficients');
    end
    
    [eqs, names, types] = equationsChemicalLog(model, state);

    species = cellfun(@(x) exp(x), logSpecies, 'UniformOutput', false);
    
    CVCind = strcmpi(model.CVC, chemsys.elementNames)';

%     %% recalculate mass balance on CVC
%     eqInd = strcmpi(names, ['Conservation of ' model.CVC]);
% 
%     masssum = 0;
%     for k = 1 : chemsys.nC
%         masssum = masssum + chemsys.compositionMatrix(CVCind,k).*species{k};
%     end
%     masssum = masssum - CVC;
% 
%     eqs{eqInd} = log(masssum) - logElements{CVCind};


    %% charge balance

    CV = chemsys.chargeVector;
    eInd = strcmpi('e-', chemsys.speciesNames);
    CV(1,eInd) = 0;

    CVp = max(CV, 0);
    CVn = max(-CV, 0);
    
    pos = 0;
    neg = 0;
    
    for k = 1 : chemsys.nC
        pos = pos + CVp(1,k).*species{k};
        neg = neg + CVn(1,k).*species{k};
    end

    eqs{end+1} = log(pos) - log(neg);
    names{end+1} = 'charge balance equation';
    types{end+1} = 'cell';
end

