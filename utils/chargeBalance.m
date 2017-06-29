function [state, model] = chargeBalance(model, state)

    comps = cell(1, model.nC);

    [comps{:}] = model.getProps(state, model.CompNames{:});
            
    T = 298;
    e_w = 87.740 - 0.4008*(T-273.15) + 9.398e-4*(T-273.15)^2 - 1.410e-6*(T-273.15)^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T)^(-3/2);
                                     

    % calculate activity
    indS = cellfun(@(x) isempty(x), regexpi(model.CompNames, '>'));
    nC = sum(indS);    
    
    model.CompActivityNames  = cellfun(@(name) ['a', name], model.CompNames(indS), ...
                                         'uniformoutput', false);
    charge = 0;
    chargeabs = 0;
    for i = 1 : nC
        charge = charge + (model.ChargeVector(1,i).*comps{i});
        chargeabs = chargeabs + abs(model.ChargeVector(1,i)).*comps{i};
    end

    state = model.setProp(state, 'chargebalance', charge./chargeabs*100);


end
