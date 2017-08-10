function [eqs, names, types] = equationsChemicalInit(logcomps, logmasterComps, comboComps, model)


    [eqs, names, types] = equationsChemicalLog(logcomps, logmasterComps, model);

    comps = cellfun(@(x) exp(x), logcomps, 'UniformOutput', false);

    %% combination matrix
    for i = 1 : model.nLC
        combSum = 0;
        for k = 1 : model.nC
            combSum = combSum + model.CombinationMatrix(i,k).*comps{k};
        end
        % log version
        eqs{end + 1} = log(combSum) - log(comboComps{i});
        % plain version
        % eqs{end + 1} = combSum - comboComps{i};
        names{end + 1} = [model.CombinationNames{i}] ;
        types{end + 1} = 'cell';
    end
            
    


end
