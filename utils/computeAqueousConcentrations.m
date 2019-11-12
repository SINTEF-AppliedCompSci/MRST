function [state, model] = computeAqueousConcentrations(model, state)

    chemsys = model.chemicalSystem;

    comps = cell(1, chemsys.nC);

    [comps{:}] = model.getProps(state, chemsys.speciesNames{:});
    
    CM = chemsys.compositionMatrix;
    
    surfInd = cellfun(@(x) ~isempty(x), regexpi(chemsys.speciesNames, '>'));
    
    CM(:,surfInd) = 0;
    
    T = 298;
    e_w = 87.740 - 0.4008*(T-273.15) + 9.398e-4*(T-273.15)^2 - 1.410e-6*(T-273.15)^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T)^(-3/2);

    % Calculate activity
    indS = cellfun(@(x) isempty(x), regexpi(chemsys.elementNames, '>'));
    nC = sum(indS);    
    
    chemsys.aqueousConcentrationNames  = cellfun(@(name) [name '(aq)'], chemsys.elementNames(indS), ...
                                         'uniformoutput', false);
    model.chemicalSystem = chemsys;
    
    totals = cell(1, nC);
                                     
    for i = 1 : nC
        totals{i} = 0;
        for j = 1 : size(CM,2)
            totals{i} = totals{i} + CM(i,j)*comps{j};
        end
    end
    % create aqueousConcentrations field, if not existing, and assign default values
    if ~isfield(state, 'aqueousConcentrations')
        elementsforsize = model.getProp(state, 'elements');
        ncells = size(elementsforsize, 1);
        clear elementsforsize
        state.aqueousConcentrations = ones(ncells, nC);
    end
    
    for i = 1 : nC
        state = model.setProp(state, chemsys.aqueousConcentrationNames{i}, totals{i});
    end

end
