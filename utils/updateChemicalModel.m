function [ state ] = updateChemicalModel(model, problem, state, state0 )

if model.nP>0
    surfParam = sum(cellfun(@(x) ~isempty(x) , regexpi(model.surfaceActivityCoefficientNames, 'psi'))); 
    if surfParam > 0
        [names, mins, maxs] = computeMaxPotential(model, state0); 
    end
end

matrixVolumeFraction = model.getProp(state,'matrixVolumeFraction');

nonLogVariables = removeLogFromNames(problem.primaryVariables); 


len = cellfun(@(x) length(x), nonLogVariables);
[~,sortInd] = sort(len(:),1, 'ascend');
pVar = nonLogVariables(sortInd);

LC = model.combinationMatrix;

for i = 1 : numel(pVar)

    p = pVar{i};
    compInd = strcmpi(p, model.componentNames);

    if any(strcmpi(p, model.masterComponentNames))
        state = model.capProperty(state, p, realmin, 2.5*mol/litre);
        
    elseif ~isempty(regexpi(p, 'psi'))
        ind = find(strcmpi(p, names));
        state = model.capProperty(state, p, mins{ind}, maxs{ind});
        
    elseif any(strcmpi(p, [model.gasNames, model.solidNames, 'fluidVolumeFraction']));
        state = model.capProperty(state, p, realmin, 1-matrixVolumeFraction);
        
    elseif any(strcmpi(p, model.combinationNames))
        ind = strcmpi(p, model.combinationNames);
        combMaxMatrix = diag(1./LC(:, ind));
        maxvals = combMaxMatrix*((state.combinationComponents)');
        maxvals = (min(maxvals))';
        state = model.capProperty(state, model.combinationNames{i}, realmin, ...
                                  maxvals);
   elseif strcmpi(p, 'CVC') 
        cvcInd = strcmpi(model.CVC, model.masterComponentNames);
        cvcVal = state.masterComponents(:,cvcInd);
        state = model.capProperty(state, p, -cvcVal*0.99, cvcVal);
        
    else
        fvf = 1 - model.getProp(state, 'matrixVolumeFraction');
        maxvals = model.maxMatrices{compInd}*((state.masterComponents)');
        maxvals = (min(maxvals))'./fvf;             
        state = model.capProperty(state, p, realmin, maxvals); 
    end

end


end

