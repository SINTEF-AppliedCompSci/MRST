function [eqs, names, types] = equationsChargeBalance(model, state)
    
    logSpecies                     = model.getProp(state, 'logSpecies');
    logElements                    = model.getProp(state, 'logElements');
    combinationComponents          = model.getProp(state, 'combinationComponents');
    logPartialPressures            = model.getProp(state, 'logPartialPressures');
    logSaturationIndicies          = model.getProp(state, 'logSaturationIndicies');
    logSurfaceActivityCoefficients = model.getProp(state, 'logSurfaceActivityCoefficients');

    
    [eqs, names, types] = equationsChemicalLog(model, state);

    species = cellfun(@(x) exp(x), logSpecies, 'UniformOutput', false);
    
    CVCind = strcmpi(model.CVC, model.elementNames)';

%     %% recalculate mass balance on CVC
%     eqInd = strcmpi(names, ['Conservation of ' model.CVC]);
% 
%     masssum = 0;
%     for k = 1 : model.nC
%         masssum = masssum + model.compositionMatrix(CVCind,k).*species{k};
%     end
%     masssum = masssum - CVC;
% 
%     eqs{eqInd} = log(masssum) - logElements{CVCind};


    %% charge balance

    CV = model.chargeVector;
    eInd = strcmpi('e-', model.speciesNames);
    CV(1,eInd) = 0;

    CVp = max(CV, 0);
    CVn = max(-CV, 0);
    
    pos = 0;
    neg = 0;
    
    for k = 1 : model.nC
        pos = pos + CVp(1,k).*species{k};
        neg = neg + CVn(1,k).*species{k};
    end

    eqs{end+1} = log(pos) - log(neg);
    names{end+1} = 'charge balance equation';
    types{end+1} = 'cell';
end

