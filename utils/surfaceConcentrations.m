function [state, model] = surfaceConcentrations(model, state)

    comps = cell(1, model.nC);

    [comps{:}] = model.getProps(state, model.componentNames{:});
    
    CM = model.compositionMatrix;
    
    surfInd = cellfun(@(x) isempty(x), regexpi(model.componentNames, '>'));
    
    CM(:,surfInd) = 0;
                                         
    indS = cellfun(@(x) isempty(x), regexpi(model.masterComponentNames, '>'));
    nC = sum(indS);    
    
    model.surfaceConcentrationNames  = cellfun(@(name) [name '_surf'], model.masterComponentNames(indS), ...
                                         'uniformoutput', false);

    totals = cell(1, nC);
                                     
    for i = 1 : nC
        totals{i} = 0;
        for j = 1 : size(CM,2)
            totals{i} = totals{i} + CM(i,j)*comps{j};
        end
    end

    for i = 1 : nC
        state = model.setProp(state, model.surfaceConcentrationNames{i}, totals{i});
    end

end
