function [eqs, names, types] = equationsCompositionGuess(comps, masterComps, comboComps, model)
    
    CNames = model.CompNames;

    compInd = cellfun(@(x) isempty(x), regexpi(CNames, 'psi'));
    CNames = CNames(compInd);
    nC = numel(CNames);

    CM = model.CompositionMatrix;
    CM(:,~compInd) = [];
            
    eqs   = cell(1, model.nMC + model.nLC);
    names = cell(1, model.nMC + model.nLC);
    types = cell(1, model.nMC + model.nLC);


    for i = 1 : model.nMC
        eqs{i} = - masterComps{i};
        for k = 1 : nC
            eqs{i} = eqs{i} + CM(i,k).*comps{k};
        end
        names{i} = ['Conservation of ', model.MasterCompNames{i}] ;
    end

    %% combination matrix
    for i = 1 : model.nLC
        j = model.nMC + i;
        combSum = 0;
        for k = 1 : model.nC
            combSum = combSum + model.CombinationMatrix(i,k).*comps{k};
        end
        eqs{j} = combSum - comboComps{i};
        names{j} = [model.CombinationNames{i}] ;
    end
    
    [types{:}] = deal('cell');
    
end

