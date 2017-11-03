function [eqs, names, types] = equationsChargeBalance_nonLog(model, state, components, masterComponents, combinationComponents,...
                 partialPressures, saturationIndicies, surfaceAcitivityCoefficients, CVC)

     logComponents = cellfun(@(x) log(x), components, 'UniformOutput', false);
     logMasterComponents = cellfun(@(x) log(x), masterComponents, 'UniformOutput', false);
     logPartialPressures = cellfun(@(x) log(x), partialPressures, 'UniformOutput', false);
     logSaturationIndicies = cellfun(@(x) log(x), saturationIndicies, 'UniformOutput', false);
     logSurfaceAcitivityCoefficients = cellfun(@(x) log(x), surfaceAcitivityCoefficients, 'UniformOutput', false);

    [eqs, names, types] = equationsChemicalLog(model, state, logComponents, logMasterComponents, combinationComponents, ...
                                                       logPartialPressures, logSaturationIndicies,logSurfaceAcitivityCoefficients);

    
    CVCind = strcmpi(model.CVC, model.elementNames)';

    %% recalculate mass balance on CVC
    eqInd = strcmpi(names, ['Conservation of ' model.CVC]);

    masssum = 0;
    for k = 1 : model.nC
        masssum = masssum + model.compositionMatrix(CVCind,k).*components{k};
    end
    masssum = masssum + CVC;

    eqs{eqInd} = log(masssum) - logMasterComponents{CVCind};


    %% charge balance

    CV = model.chargeVector;
    eInd = strcmpi('e-', model.speciesNames);
    CV(1,eInd) = 0;

    CVp = max(CV, 0);
    CVn = max(-CV, 0);
    
    pos = 0;
    neg = 0;
    
    for k = 1 : model.nC
        pos = pos + CVp(1,k).*components{k};
        neg = neg + CVn(1,k).*components{k};
    end

    eqs{end+1} = log(pos) - log(neg);
    names{end+1} = 'charge balance equation';
    types{end+1} = 'cell';
end

