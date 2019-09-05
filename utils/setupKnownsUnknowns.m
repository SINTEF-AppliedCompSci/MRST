function [ model ] = setupKnownsUnknowns( model )
    % setup unknownNames
    unknownNames = horzcat(model.componentNames, model.masterComponentNames, model.combinationNames, model.solidNames, model.gasNames, model.surfaceActivityCoefficientNames, 'fluidVolumeFraction');
    ind = cellfun(@(name)(strcmpi(name, model.inputNames)), unknownNames, ...
                  'Uniformoutput', false);

    ind = cell2mat(ind');
    ind = sum(ind, 2);
    ind = logical(ind);
    model.unknownNames = unknownNames(~ind);


end

