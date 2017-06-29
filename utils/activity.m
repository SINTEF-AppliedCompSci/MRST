function [state, model] = activity(model, state)

    comps = cell(1, model.nC);

    [comps{:}] = model.getProps(state, model.CompNames{:});
            
    try 
        T = model.getProp(state, 'temperature');
    catch
        T = 298;
    end
    e_w = 87.740 - 0.4008*(T-273.15) + 9.398e-4*(T-273.15).^2 - 1.410e-6*(T-273.15).^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T).^(-3/2);
                                     

    % calculate activity
    ionDum = 0;
    indS = cellfun(@(x) isempty(x), regexpi(model.CompNames, '>'));
    nC = sum(indS);    

    model.CompActivityNames  = cellfun(@(name) ['a', name], model.CompNames(indS), ...
                                         'uniformoutput', false);
                                     
    for i = 1 : nC
        ionDum = ionDum + (model.ChargeVector(1,i).^2.*comps{i}).*litre/mol;
    end
    ion = cell(1,model.nC);
    [ion{:}] = deal((1/2)*ionDum);
    
    pg = cell(1,nC);
    for i = 1 : nC
        pg{i} = log(10).*-A.*model.ChargeVector(1,i)'.^2 .* (ion{i}.^(1/2)./(1 + ion{i}.^(1/2)) - 0.3.*ion{i});
%         if isfield(pg{i}, 'val')
%             doub = pg{i}.val;
%         else
%             doub = pg{i};
%         end
%         
%         if sum(doub) == 0 && ~strcmpi(model.CompNames{i}, 'H2O')
%             pg{i} = log(10^0.010)*ones(size(pg{i},1),1);
%         end
    end

    % Reaction matrix, activities only apply to laws of mass action

    for k = 1 : nC
        activities{k} =  (exp(pg{k}) .* comps{k});
        state = model.setProp(state, model.CompActivityNames{k}, activities{k});
    end



end
