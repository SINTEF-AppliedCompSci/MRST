function [eqs, names, types] = equationsChemicalInit(model, state, logFluidVolumeFraction, logComponents, logMasterComponents, combinationComponents, ...
                                                       logGasVolumeFraction, logSolidVolumeFraction,logSurfaceAcitivityCoefficients)


    [eqs, names, types] = equationsChemicalLog(model, state, logFluidVolumeFraction, logComponents, logMasterComponents, combinationComponents, ...
                                                       logGasVolumeFraction, logSolidVolumeFraction,logSurfaceAcitivityCoefficients);

    components = cellfun(@(x) exp(x), logComponents, 'UniformOutput', false);

    %% combination matrix
    for i = 1 : model.nLC
        combSum = 0;
        for k = 1 : model.nC
            combSum = combSum + model.combinationMatrix(i,k).*components{k};
        end
        eqs{end + 1} = log(combSum) - log(combinationComponents{i});
        names{end + 1} = [model.combinationNames{i}] ;
        types{end + 1} = 'cell';
    end
            
    


end
