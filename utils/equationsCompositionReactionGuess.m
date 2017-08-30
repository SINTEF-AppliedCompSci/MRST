function [eqs, names, types] = equationsCompositionReactionGuess(logcomps, logmasterComps, comboComps, model)
    
%             comps = cellfun(@(x) x*litre/mol, comps,'UniformOutput', false);
%             masterComps = cellfun(@(x) x*litre/mol, masterComps,'UniformOutput', false);
%             
% begin solving equation

    CM = model.CompositionMatrix;
    RM = model.ReactionMatrix;

    nC = size(CM,2);
    
    comps = cellfun(@(x) exp(x), logcomps, 'UniformOutput', false);
    
    logK = model.LogReactionConstants;
    
    eqs   = cell(1, model.nMC + model.nR + model.nLC);
    names = cell(1, model.nMC + model.nR + model.nLC);
    types = cell(1, model.nMC + model.nR + model.nLC);

    %% reaction matrix
    for i = 1 : model.nR


        eqs{i} = -logK(i);
        for k = 1 : nC
            eqs{i} = eqs{i} + RM(i, k).*logcomps{k};
        end
        names{i} = model.rxns{i};
    end
    
    %% composition matrix
    for i = 1 : model.nMC
        j = model.nR + i;
        masssum = 0;
        
        for k = 1 : nC
            masssum = masssum + CM(i,k).*comps{k};
        end
        
        eqs{j} = log(masssum) - logmasterComps{i};
        
        names{j} = ['Conservation of ', model.MasterCompNames{i}] ;
    end

    %% combination matrix
    for i = 1 : model.nLC
        j = model.nR + model.nMC + i;
        combSum = 0;
        for k = 1 : model.nC
            combSum = combSum + model.CombinationMatrix(i,k).*comps{k};
        end
        eqs{j} = combSum - comboComps{i};
        names{j} = [model.CombinationNames{i}] ;
    end
            
    [types{:}] = deal('cell');
    
end

