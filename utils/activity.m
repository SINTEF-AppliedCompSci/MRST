function [state, model] = activity(model, state)

    components = cell(1, model.nC);

    [components{:}] = model.getProps(state, model.speciesNames{:});
            
    T = model.getProp(state, 'temperature');
    e_w = 87.740 - 0.4008*(T-273.15) + 9.398e-4*(T-273.15).^2 - 1.410e-6*(T-273.15).^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T).^(-3/2);
                                     

    % calculate activity
    ionDum = 0;
    indS = cellfun(@(x) isempty(x), regexpi(model.speciesNames, '>'));
    nC = sum(indS);    

    model.activityNames  = cellfun(@(name) ['a', name], model.speciesNames(indS), ...
                                         'uniformoutput', false);
    
    CV = model.chargeVector;
    eInd = strcmpi('e-', model.speciesNames);
    CV(1,eInd) = 0;
    
    for i = 1 : nC
        ionDum = ionDum + (CV(1,i).^2.*components{i}).*litre/mol;
    end
    ion = cell(1,model.nC);
    [ion{:}] = deal((1/2)*ionDum);
    
    pg = cell(1,nC);
    for i = 1 : nC
        pg{i} = log(10).*-A.*CV(1,i)'.^2 .* (ion{i}.^(1/2)./(1 + ion{i}.^(1/2)) - 0.3.*ion{i});
        if CV(1,i) == 0
            pg{i} = ion{i}*0.1;
        end
    end

    % Reaction matrix, activities only apply to laws of mass action

    for k = 1 : nC
        activities{k} =  (exp(pg{k}) .* components{k});
        state = model.setProp(state, model.activityNames{k}, activities{k});
    end



end
