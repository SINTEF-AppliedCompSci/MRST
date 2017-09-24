function [eqs, names, types] = equationsCompositionReactionGuess(model, state, logFluidVolumeFraction, logComponents, logMasterComponents, combinationComponents, logGasVolumeFractions, logSolidVolumeFractions)
    
%     if model.nG > 0
%         partialPressures = cell(1,model.nG);
%         [partialPressures{:}] = deal(model.getProps(state, model.partialPressureNames{:}));
%         logPartialPressures = cellfun(@(x) log(x), partialPressures, 'UniformOutput', false);
%     end

    if model.nS > 0
        solidDensities = cell(1,model.nS);
        [solidDensities{:}] = model.getProps(state,  model.solidDensityNames{:});
    end
    
    matrixVolumeFraction = model.getProps(state,  'matrixVolumeFraction');
    T = model.getProps(state, 'temperature');
    R   = 8.3144621;             	% Gas Constant [J/(K mol)]

        
    CM = model.compositionMatrix;
    RM = model.reactionMatrix;

    nC = size(CM,2);
    
    components       = cellfun(@(x) exp(x), logComponents, 'UniformOutput', false);
    gasVolumeFractions    = cellfun(@(x) exp(x), logGasVolumeFractions, 'UniformOutput', false);
    solidVolumeFractions  = cellfun(@(x) exp(x), logSolidVolumeFractions, 'UniformOutput', false);
    masterComponents = cellfun(@(x) exp(x), logMasterComponents, 'UniformOutput', false);
    
    fluidVolumeFraction = exp(logFluidVolumeFraction);
    
%     gasVolumes = cellfun(@(x) exp(x), logGasVolumes, 'UniformOutput', false);
    
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
%         
%         for k = 1 : model.nG
%             eqs{i} = eqs{i} + model.gasReactionMatrix(i,k).*logGasComps{k};
%         end
        
        names{i} = model.rxns{i};
    end
    
    %% composition matrix
    for i = 1 : model.nMC
        
        j = model.nR + i;
        masssum = 0;
        
        for k = 1 : nC
            masssum = masssum + CM(i,k).*components{k}.*fluidVolumeFraction;
        end
%         
%         for k = 1 : model.nG
%             masssum = masssum + model.gasContributionMatrix(i,k).*gasVolumeFractions{k}.*partialPressures{k}./(R*T);
%         end
        
        for k = 1 : model.nS
            masssum = masssum + model.solidContributionMatrix(i,k).*solidVolumeFractions{k}.*solidDensities{k};
        end
        
        eqs{j} = logMasterComponents{i} - log(masssum);
                
        names{j} = ['Conservation of ', model.masterComponentNames{i}] ;
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
    
    %% conservation of volume
    
    vol = fluidVolumeFraction + matrixVolumeFraction;
    for i = 1 : model.nS
        vol = vol + solidVolumeFractions{i};
    end
%     for i = 1 : model.nG
%         vol = vol + gasVolumeFractions{i};
%     end    
    
    eqs{end+1} = log(1) - log(vol);
    names{end+1} = 'Conservation of volume';
    types{end+1} = [];
    
    [types{:}] = deal('cell');
    
    
    
end

