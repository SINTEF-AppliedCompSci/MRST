function [eqs, names, types] = equationsChemicalInit(model, state, logFluidVolumeFraction, logSpecies, logElements, combinationComponents, ...
                                                       logGasVolumeFraction, logSolidVolumeFraction,logSurfaceActivityCoefficients)


    [eqs, names, types] = equationsChemicalLog(model, state, logFluidVolumeFraction, logSpecies, logElements, combinationComponents, ...
                                                       logGasVolumeFraction, logSolidVolumeFraction,logSurfaceActivityCoefficients);

    species = cellfun(@(x) exp(x), logSpecies, 'UniformOutput', false);

    %% combination matrix
    for i = 1 : chemsys.nLC
        combSum = 0;
        for k = 1 : chemsys.nC
            combSum = combSum + chemsys.combinationMatrix(i,k).*species{k};
        end
        eqs{end + 1} = log(combSum) - log(combinationComponents{i});
        names{end + 1} = [chemsys.combinationNames{i}] ;
        types{end + 1} = 'cell';
    end
            
    


end
