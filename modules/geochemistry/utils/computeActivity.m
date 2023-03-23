function [state, model] = computeActivity(model, state)

    chemsys = model.chemicalSystem;
    
    species = cell(1, chemsys.nC);

    [species{:}] = model.getProps(state, chemsys.speciesNames{:});
            
    T = model.getProp(state, 'temperature');
    e_w = 87.740 - 0.4008*(T-273.15) + 9.398e-4*(T-273.15).^2 - 1.410e-6*(T-273.15).^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T).^(-3/2);
                                     

    % calculate activity
    ionDum = 0;
    indS = cellfun(@(x) isempty(x), regexpi(chemsys.speciesNames, '>'));
    nC = sum(indS);    

    chemsys.activityNames  = cellfun(@(name) ['a', name], chemsys.speciesNames(indS), ...
                                   'uniformoutput', false);
    model.chemicalSystem = chemsys;
    
    CV = chemsys.chargeVector;
    eInd = strcmpi('e-', chemsys.speciesNames);
    CV(1,eInd) = 0;
    
    for i = 1 : nC
        ionDum = ionDum + (CV(1,i).^2.*species{i}).*litre/mol;
    end
    ion = cell(1,chemsys.nC);
    [ion{:}] = deal((1/2)*ionDum);
    
    pg = cell(1, nC);
    for i = 1 : nC
        pg{i} = log(10).*-A.*CV(1,i)'.^2 .* (ion{i}.^(1/2)./(1 + ion{i}.^(1/2)) - 0.3.*ion{i});
        if CV(1,i) == 0
            pg{i} = ion{i}*0.1;
        end
    end

    % we check if activity field exits, if not we create it
    if ~isfield(state, 'activities')
        % get number of rows (number of evaluated cells) from element
        elementsforsize = model.getProp(state, 'elements');
        ncells = size(elementsforsize, 1);
        clear elementsforsize
        ncomp = numel(chemsys.activityNames);
        state.activities = ones(ncells, ncomp);
    end
        
    % Reaction matrix, activities only apply to laws of mass action
    for k = 1 : nC
        activity = (exp(pg{k}) .* species{k});
        state = model.setProp(state, chemsys.activityNames{k}, activity);
    end

end
