function [eqs, names, types] = equationsChargeBalance(model, state, logComponents, logMasterComponents, combinationComponents,...
                 logPartialPressures, logSaturationIndicies, logSurfaceAcitivityCoefficients)

    [eqs, names, types] = equationsChemicalLog(model, state, logComponents, logMasterComponents, combinationComponents, ...
                                                       logPartialPressures, logSaturationIndicies,logSurfaceAcitivityCoefficients);

    components = cellfun(@(x) exp(x), logComponents, 'UniformOutput', false);
    
    CVCind = strcmpi(model.CVC, model.elementNames)';

%     %% recalculate mass balance on CVC
%     eqInd = strcmpi(names, ['Conservation of ' model.CVC]);
% 
%     masssum = 0;
%     for k = 1 : model.nC
%         masssum = masssum + model.compositionMatrix(CVCind,k).*components{k};
%     end
%     masssum = masssum - CVC;
% 
%     eqs{eqInd} = log(masssum) - logMasterComponents{CVCind};


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

