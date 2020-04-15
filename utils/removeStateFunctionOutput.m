function model = removeStateFunctionOutput(model)
    if isprop(model, 'parentModel')
        model.parentModel = removeStateFunctionOutput(model.parentModel);
    end
    model.OutputStateFunctions = {};
end