function [ state ] = updateChemicalModel(model, problem, state, state0 )

if model.nP>0
    surfParam = sum(cellfun(@(x) ~isempty(x) , regexpi(model.surfaceActivityCoefficientNames, 'psi'))); 
    if surfParam > 0
        [names, mins, maxs] = computeMaxPotential(model, state0); 
    end
end

nonLogVariables = removeLogFromNames(problem.primaryVariables); 


len = cellfun(@(x) length(x), nonLogVariables);
[~,sortInd] = sort(len(:),1, 'ascend');
pVar = nonLogVariables(sortInd);

LC = model.combinationMatrix;
CM = model.compositionMatrix;

if isfield(state, 'CVC')
    mult = 100;
end

for i = 1 : numel(pVar)

    p = pVar{i};
    compInd = strcmpi(p, model.speciesNames);

    if any(strcmpi(p, model.elementNames))
        state = model.capProperty(state, p, realmin, 300*mol/litre);
        
    elseif ~isempty(regexpi(p, 'psi'))
        ind = find(strcmpi(p, names));
        state = model.capProperty(state, p, mins{ind}, maxs{ind});
        
    elseif any(strcmpi(p, model.combinationNames))
        ind = strcmpi(p, model.combinationNames);
        
        conMat = CM;
        conMat(:,LC(ind,:)==0) = 0; 
        conMat = sum(conMat,2);
        maxvals = state.elements*conMat;
        state = model.capProperty(state, p, -maxvals, ...
                                  maxvals);
                              
%    elseif strcmpi(p, 'CVC') 
%         cvcInd = strcmpi(model.CVC, model.elementNames);
%         cvcVal = state.elements(:,cvcInd);
%         figure(1); hold on;
%         plot(model.getProp(state, p),'-k')
%         state = model.capProperty(state, p, -0.99*cvcVal, mult*cvcVal);
%         plot(model.getProp(state, p),'--r')
%         if true
%         end
        
    elseif any(strcmpi(p, model.speciesNames))
        maxvals = model.maxMatrices{compInd}*((state.elements)');
        maxvals = (min(maxvals))';
%         if ~isempty(strfind(p,model.CVC)) 
%             state = model.capProperty(state, p, 1e-100, mult*maxvals);
%         else 
            state = model.capProperty(state, p, 1e-100, maxvals);
%         end
    end

end

if true
end

end

