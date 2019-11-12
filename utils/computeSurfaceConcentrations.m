function [state, model] = computeSurfaceConcentrations(model, state)

    chemsys = model.chemicalSystem;

    comps = cell(1, chemsys.nC);

    [comps{:}] = model.getProps(state, chemsys.speciesNames{:});
    
    CM = chemsys.compositionMatrix;
    
    surfInd = cellfun(@(x) isempty(x), regexpi(chemsys.speciesNames, '>'));
    
    CM(:,surfInd) = 0;
                                         
    indS = cellfun(@(x) isempty(x), regexpi(chemsys.elementNames, '>'));
    nC = sum(indS);    
    
    chemsys.surfaceConcentrationNames  = cellfun(@(name) [name '(surf)'], chemsys.elementNames(indS), ...
                                         'uniformoutput', false);
    model.chemicalSystem = chemsys;

    totals = cell(1, nC);
                                     
    for i = 1 : nC
        totals{i} = 0;
        for j = 1 : size(CM,2)
            totals{i} = totals{i} + CM(i,j)*comps{j};
        end
    end

    % create surfaceConcentrations field, if not existing, and assign default values
    if ~isfield(state, 'surfaceConcentrations')
        elementsforsize = model.getProp(state, 'elements');
        ncells = size(elementsforsize, 1);
        clear elementsforsize;
        state.surfaceConcentrations = ones(ncells, nC);
    end
    
    for i = 1 : nC
        state = model.setProp(state, chemsys.surfaceConcentrationNames{i}, totals{i});
    end

end
