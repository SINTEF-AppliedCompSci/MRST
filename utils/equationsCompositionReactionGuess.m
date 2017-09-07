function [eqs, names, types] = equationsCompositionReactionGuess(state, logporo, logcomps, logmasterComps, comboComps, logGasComps, logSolidComps, model)
    
%     if model.nG > 0
%         partialPressures = cell(1,model.nG);
%         [partialPressures{:}] = deal(model.getProps(state, model.partialPressureNames{:}));
%         logPartialPressures = cellfun(@(x) log(x), partialPressures, 'UniformOutput', false);
%     end

    if model.nS > 0
        solidDensities = cell(1,model.nS);
        [solidDensities{:}] = deal(model.getProps(state,  model.solidDensityNames{:}));
    end
    
    T = model.getProps(state, 'temp');
    
    R   = 8.3144621;             	% Gas Constant [J/(K mol)]

        
    CM = model.CompositionMatrix;
    RM = model.ReactionMatrix;

    nC = size(CM,2);
    
    comps = cellfun(@(x) exp(x), logcomps, 'UniformOutput', false);
    gasComps = cellfun(@(x) exp(x), logGasComps, 'UniformOutput', false);
    solidComps = cellfun(@(x) exp(x), logSolidComps, 'UniformOutput', false);
    poro = exp(logporo);
    
%     gasVolumes = cellfun(@(x) exp(x), logGasVolumes, 'UniformOutput', false);
    
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
        
        for k = 1 : model.nG
            eqs{i} = eqs{i} + model.GasReactionMatrix(i,k).*logGasComps{k};
        end
        
        names{i} = model.rxns{i};
    end
    
    %% composition matrix
    for i = 1 : model.nMC
        
        j = model.nR + i;
        masssum = 0;
        
        for k = 1 : nC
            masssum = masssum + CM(i,k).*comps{k}.*poro;
        end
%         
%         for k = 1 : model.nG
%             masssum = masssum + model.GasCompMatrix(i,k).*gasComps{k}.*partialPressures{k}./(R*T);
%         end
        
        for k = 1 : model.nS
            masssum = masssum + model.SolidCompMatrix(i,k).*solidComps{k}.*solidDensities{k};
        end
        
        eqs{j} = log(masssum) - (logmasterComps{i} + logporo);
                
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
    
    %% conservation of volume
    
    vol = poro;
    for i = 1 : model.nS
        vol = vol + solidComps{i};
    end
%     for i = 1 : model.nG
%         vol = vol + gasComps{i};
%     end    
    
    eqs{end+1} = log(1) - log(vol);
    names{end+1} = 'Conservation of volume';
    types{end+1} = [];
    
    [types{:}] = deal('cell');
    
    
    
end

