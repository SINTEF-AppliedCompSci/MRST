function [eqs, names, types] = equationsCompositionReactionGuess(model, state, logComponents, logMasterComponents, combinationComponents, logPartialPressures, logSaturationIndicies)
    
   
    CM = model.compositionMatrix;
    RM = model.reactionMatrix;
    GM =  model.gasReactionMatrix;
    SM = model.solidReactionMatrix;
    
    nC = size(CM,2);
    
    components       = cellfun(@(x) exp(x), logComponents, 'UniformOutput', false);
        
    logK = model.logReactionConstants;
    
    eqs   = cell(1, model.nMC + model.nR + model.nLC);
    names = cell(1, model.nMC + model.nR + model.nLC);
    types = cell(1, model.nMC + model.nR + model.nLC);

    %% reaction matrix
    for i = 1 : model.nR

        eqs{i} = -logK(i);
        
        for k = 1 : nC
            eqs{i} = eqs{i} + RM(i, k).*logComponents{k};
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
        eqs{j} = combSum - combinationComponents{i};
        names{j} = [model.combinationNames{i}] ;
    end
    
    [types{:}] = deal('cell');
    
    
    
end

