function [state, model] = computeChargeBalance(model, state)

    chemsys = model.chemicalSystem;
    
    species = cell(1, chemsys.nC);

    [species{:}] = model.getProps(state, chemsys.speciesNames{:});
           
    CV = chemsys.chargeVector;
    eInd = strcmpi('e-', chemsys.speciesNames);
    CV(1,eInd) = 0;

    charge = 0;
    chargeabs = 0;
    for i = 1 : chemsys.nC
        charge = charge + CV(1,i).*species{i};
        chargeabs = chargeabs + abs(CV(1,i)).*species{i};
    end

    state = model.setProp(state, 'chargebalance', charge./chargeabs);


end
