function [state, primaryVars] = setupDynamicStateOilWater(model, state, useAD, ...
                                                                 reverseMode)
    % Properties at current timestep
    [p, sW, wellSol] = model.getProps(state, ...
        'pressure', 'water', 'wellSol');

    [wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);


    primaryVars = {};
    if useAD
        [p, sW, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, sW, wellVars{:});
        if reverseMode
            primaryVars = {'pressure', 'sW'};
        else
            primaryVars = {'pressure', 'sW', wellVarNames{:}};
        end
    end
    sO = 1 - sW;

    sat = {sW, sO};
    wellSol = DynamicState(wellSol, [wellVarNames, 'wellmap'], [wellVars, wellMap]);
    
    [fp, fpname] = model.FlowPropertyFunctions.getPropertyContainer();
    
    state = DynamicState(state, {'pressure', 's', 'wellSol', fpname},...
                                {p, sat, wellSol, fp});
end