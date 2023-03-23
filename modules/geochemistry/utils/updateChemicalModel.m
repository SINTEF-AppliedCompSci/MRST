function state  = updateChemicalModel(model, problem, state, state0 )

    chemsys = model.chemicalSystem;
    
    if chemsys.nP>0
        surfParam = sum(cellfun(@(x) ~isempty(x) , regexpi(chemsys.surfaceActivityCoefficientNames, 'psi'))); 
        if surfParam > 0
            [names, mins, maxs] = computeMaxPotential(model, state0); 
        end
    end

    nonLogVariables = removeLogFromNames(problem.primaryVariables); 

    len = cellfun(@(x) length(x), nonLogVariables);
    [~,sortInd] = sort(len(:),1, 'ascend');
    pVar = nonLogVariables(sortInd);

    LC = chemsys.combinationMatrix;
    CM = chemsys.compositionMatrix;


    for i = 1 : numel(pVar)

        p = pVar{i};
        compInd = strcmpi(p, chemsys.speciesNames);

        if any(strcmpi(p, chemsys.elementNames))
            state = model.capProperty(state, p, realmin, 300*mol/litre);
            
        elseif ~isempty(regexpi(p, 'psi'))
            ind = find(strcmpi(p, names));
            state = model.capProperty(state, p, mins{ind}, maxs{ind});
            
        elseif any(strcmpi(p, chemsys.combinationNames))
            ind = strcmpi(p, chemsys.combinationNames);
            
            conMat = CM;
            conMat(:,LC(ind,:)==0) = 0; 
            conMat = sum(conMat,2);
            maxvals = state.elements*conMat;
            state = model.capProperty(state, p, -maxvals, maxvals);
            
        elseif any(strcmpi(p, chemsys.speciesNames))
            maxvals = chemsys.maxMatrices{compInd}*((state.elements)');
            maxvals = (min(maxvals))';
            state = model.capProperty(state, p, realmin, maxvals);
        end

    end


end

