function [state, model] = aqueousConcentrations(model, state)

    comps = cell(1, model.nC);

    [comps{:}] = model.getProps(state, model.componentNames{:});
    
    CM = model.compositionMatrix;
    
    surfInd = cellfun(@(x) ~isempty(x), regexpi(model.componentNames, '>'));
    
    CM(:,surfInd) = 0;
    
    T = 298;
    e_w = 87.740 - 0.4008*(T-273.15) + 9.398e-4*(T-273.15)^2 - 1.410e-6*(T-273.15)^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T)^(-3/2);
                                     

    % calculate activity
    indS = cellfun(@(x) isempty(x), regexpi(model.masterComponentNames, '>'));
    nC = sum(indS);    
    
    model.aqueousConcentrationNames  = cellfun(@(name) [name '_aq'], model.masterComponentNames(indS), ...
                                         'uniformoutput', false);

    totals = cell(1, nC);
                                     
    for i = 1 : nC
        totals{i} = 0;
        for j = 1 : size(CM,2)
            totals{i} = totals{i} + CM(i,j)*comps{j};
        end
    end

    for i = 1 : nC
        state = model.setProp(state, model.aqueousConcentrationNames{i}, totals{i});
    end

end
