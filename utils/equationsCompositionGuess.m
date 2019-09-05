function [eqs, names, types] = equationsCompositionGuess(poro, comps, masterComps, comboComps, gasComps, solidComps, model)
    
    CNames = model.CompNames;

    compInd = cellfun(@(x) isempty(x), regexpi(CNames, 'psi'));
    CNames = CNames(compInd);
    nC = numel(CNames);

    CM = model.CompositionMatrix;
    CM(:,~compInd) = [];
            
    eqs   = cell(1, model.nMC + model.nLC);
    names = cell(1, model.nMC + model.nLC);
    types = cell(1, model.nMC + model.nLC);


    %% composition matrix
    for i = 1 : model.nMC
        eqs{i} = - masterComps{i}.*poro;
        for k = 1 : nC
            eqs{i} = eqs{i} + CM(i,k).*comps{k}.*poro;
        end
        
        % gasComps and solidComps have units of volume
        for k = 1 : model.nG
            eqs{i} = eqs{i} + model.GasCompMatrix(i,k).*gasComps{k};
        end
        for k = 1 : model.nS
            eqs{i} = eqs{i} + model.SolidCompMatrix(i,k).*solidComps{k}.*model.solidDensities(k);
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
    
    %% conservation of volume
    vol = 1 - poro;
    for i = 1 : model.nS
        vol = vol - solidComps{i};
    end
    for i = 1 : model.nG
        vol = vol - gasComps{i};
    end    
    
    eqs{end+1} = vol;
    names{end+1} = 'Conservation of volume';
    types{end+1} = [];
    
    %%
    
    [types{:}] = deal('cell');
    
end

