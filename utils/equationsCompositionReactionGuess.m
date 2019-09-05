function [eqs, names, types] = equationsCompositionReactionGuess(model, state, logComponents, logMasterComponents, combinationComponents, logPartialPressures, logSaturationIndicies)
    
   
    CM = model.compositionMatrix;
    RM = model.reactionMatrix;
    GM =  model.gasReactionMatrix;
    SM = model.solidReactionMatrix;
    
    nC = size(CM,2);
    
    components       = cellfun(@(x) exp(x), logComponents, 'UniformOutput', false);
    masterComponents = cellfun(@(x) exp(x), logMasterComponents,'UniformOutput', false);

    logK = model.logReactionConstants;
    
    eqs   = cell(1, model.nMC + model.nR + model.nLC);
    names = cell(1, model.nMC + model.nR + model.nLC);
    types = cell(1, model.nMC + model.nR + model.nLC);

    %% mol fractions
    surfMat = repmat(model.surfMaster, 1, model.nC).*CM;
    surfTest = logical(sum(surfMat));
    
    moleFraction = components;
    
    for i = 1 : model.nC
        if surfTest(i)
            surfDen = 0;
            surfNum = 0;
            for j = 1 : model.nMC
                surfNum = surfNum + CM(j,i).*model.surfMaster(j);
                surfDen = surfDen + double(logical(CM(j,i).*model.surfMaster(j)))*masterComponents{j};
            end
            moleFraction{i} = (surfNum./surfDen).*components{i};
        end

    end
    logMoleFraction = cellfun(@(x) log(x), moleFraction, 'UniformOutput', false);
    
    %% reaction matrix
    for i = 1 : model.nR

        eqs{i} = -logK{i}(:);
        
        for k = 1 : nC
            eqs{i} = eqs{i} + RM(i, k).*logMoleFraction{k};
        end
        
        % gas reactions
        for k = 1 : model.nG
            eqs{i} = eqs{i} + GM(i,k).*logPartialPressures{k};
        end
        
        % solid reactions
        for k = 1 : model.nS
            eqs{i} = eqs{i} + SM(i,k).*logSaturationIndicies{k};
        end
        
        names{i} = model.rxns{i};
    end
    
    %% composition matrix
    for i = 1 : model.nMC
        
        j = model.nR + i;
        masssum = 0;
        
        for k = 1 : nC
            masssum = masssum + CM(i,k).*components{k};
        end

        eqs{j} = logMasterComponents{i} - log(masssum);
                
        names{j} = ['Conservation of ', model.elementNames{i}] ;
    end

    %% combination matrix
    for i = 1 : model.nLC
        j = model.nR + model.nMC + i;
        combSum = 0;
        for k = 1 : model.nC
            combSum = combSum + model.combinationMatrix(i,k).*components{k};
        end
        if any(combSum<=0) || any(combinationComponents{i}<=0)
            eqs{j} = combSum - combinationComponents{i};
        else
            eqs{j} = log(combSum) - log(combinationComponents{i});
        end
        
        names{j} = [model.combinationNames{i}] ;
    end
    
    [types{:}] = deal('cell');
    
    
    
end

