function [state, model] = chargeBalance(model, state)

    components = cell(1, model.nC);

    [components{:}] = model.getProps(state, model.speciesNames{:});
           
    CV = model.chargeVector;
    eInd = strcmpi('e-', model.speciesNames);
    CV(1,eInd) = 0;

    charge = 0;
    chargeabs = 0;
    for i = 1 : model.nC
        charge = charge + CV(1,i).*components{i};
        chargeabs = chargeabs + abs(CV(1,i)).*components{i};
    end

    state = model.setProp(state, 'chargebalance', charge./chargeabs);


end
