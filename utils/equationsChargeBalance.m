function [eqs, names, types] = equationsChargeBalance(model, state, logComponents, logMasterComponents, combinationComponents,...
                 logGasVolumeFractions, logSolidVolumeFractions, logFluidVolumeFraction, logSurfaceAcitivityCoefficients, CVC)
    

    fluidVolumeFraction = exp(logFluidVolumeFraction);
     components = cellfun(@(x) exp(x), logComponents, 'UniformOutput', false);
%     logMasterComponents = cellfun(@(x) log(x), masterComponents, 'UniformOutput', false);
%     logSolidVolumeFractions = cellfun(@(x) log(x), solidVolumeFractions, 'UniformOutput', false);
%     logGasVolumeFractions = cellfun(@(x) log(x), gasVolumeFractions, 'UniformOutput', false);
%     logSurfaceAcitivityCoefficients = cellfun(@(x) log(x), surfaceAcitivityCoefficients, 'UniformOutput', false);

    [eqs, names, types] = equationsChemicalLog(model, state, logFluidVolumeFraction, logComponents, logMasterComponents, combinationComponents, ...
                                                       logGasVolumeFractions, logSolidVolumeFractions,logSurfaceAcitivityCoefficients);


    CVCind = strcmpi(model.CVC, model.masterComponentNames)';

    %% recalculate mass balance on CVC
    eqInd = strcmpi(names, ['Conservation of ' model.CVC]);

    masssum = 0;
    for k = 1 : model.nC
        masssum = masssum + model.compositionMatrix(CVCind,k).*components{k}.*fluidVolumeFraction;
    end
    masssum = masssum+CVC;

    eqs{eqInd} = log(masssum) - logMasterComponents{CVCind};


    %% charge balance

    CV = model.chargeVector;
    eInd = strcmpi('e-', model.componentNames);
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

