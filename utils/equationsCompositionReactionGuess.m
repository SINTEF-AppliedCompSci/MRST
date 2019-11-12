function [eqs, names, types] = equationsCompositionReactionGuess(model, state)
    
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
    
    CM = chemsys.compositionMatrix;
    RM = chemsys.reactionMatrix;
    GM = chemsys.gasReactionMatrix;
    SM = chemsys.solidReactionMatrix;
    
    nC = size(CM,2);
    
    species  = cellfun(@(x) exp(x), logSpecies, 'UniformOutput', false);
    elements = cellfun(@(x) exp(x), logElements,'UniformOutput', false);

    logK = chemsys.logReactionConstants;
    
    eqs   = cell(1, chemsys.nMC + chemsys.nR + chemsys.nLC);
    names = cell(1, chemsys.nMC + chemsys.nR + chemsys.nLC);
    types = cell(1, chemsys.nMC + chemsys.nR + chemsys.nLC);

    %% mol fractions
    surfMat = repmat(chemsys.surfMaster, 1, chemsys.nC).*CM;
    surfTest = logical(sum(surfMat));
    
    moleFraction = species;
    
    for i = 1 : chemsys.nC
        if surfTest(i)
            surfDen = 0;
            surfNum = 0;
            for j = 1 : chemsys.nMC
                surfNum = surfNum + CM(j,i).*chemsys.surfMaster(j);
                surfDen = surfDen + double(logical(CM(j,i).*chemsys.surfMaster(j)))*elements{j};
            end
            moleFraction{i} = (surfNum./surfDen).*species{i};
        end

    end
    logMoleFraction = cellfun(@(x) log(x), moleFraction, 'UniformOutput', false);
    
    %% reaction matrix
    for i = 1 : chemsys.nR

        eqs{i} = -logK{i}(:);
        
        for k = 1 : nC
            eqs{i} = eqs{i} + RM(i, k).*logMoleFraction{k};
        end
        
        % gas reactions
        for k = 1 : chemsys.nG
            eqs{i} = eqs{i} + GM(i,k).*logPartialPressures{k};
        end
        
        % solid reactions
        for k = 1 : chemsys.nS
            eqs{i} = eqs{i} + SM(i,k).*logSaturationIndicies{k};
        end
        
        names{i} = chemsys.rxns{i};
    end
    
    %% composition matrix
    for i = 1 : chemsys.nMC
        
        j = chemsys.nR + i;
        masssum = 0;
        
        for k = 1 : nC
            masssum = masssum + CM(i,k).*species{k};
        end

        eqs{j} = logElements{i} - log(masssum);
                
        names{j} = ['Conservation of ', chemsys.elementNames{i}] ;
    end

    %% combination matrix
    for i = 1 : chemsys.nLC
        j = chemsys.nR + chemsys.nMC + i;
        combSum = 0;
        for k = 1 : chemsys.nC
            combSum = combSum + chemsys.combinationMatrix(i,k).*species{k};
        end
        if any(combSum<=0) || any(combinationComponents{i}<=0)
            eqs{j} = combSum - combinationComponents{i};
        else
            eqs{j} = log(combSum) - log(combinationComponents{i});
        end
        
        names{j} = [chemsys.combinationNames{i}] ;
    end
    
    [types{:}] = deal('cell');
    
    
    
end

