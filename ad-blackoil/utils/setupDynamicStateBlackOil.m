function [state, primaryVars] = setupDynamicStateBlackOil(model, state, useAD)
    % Properties at current timestep
    [p, sW, sG, rs, rv, wellSol] = model.getProps(state, ...
        'pressure', 'water', 'gas', 'rs', 'rv', 'wellSol');

    [wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);


    %Initialization of primary variables ----------------------------------
    st  = model.getCellStatusVO(state,  1-sW-sG, sW, sG);
    if model.disgas || model.vapoil
        % X is either Rs, Rv or Sg, depending on each cell's saturation status
        x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
        gvar = 'x';
    else
        x = sG;
        gvar = 'sG';
    end

    if useAD
        % define primary varible x and initialize
        if model.water
            [p, sW, x, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, sW, x, wellVars{:});
        else
            [p, x, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, x, wellVars{:});
        end
    end
    if model.water
        s_hydrocarbon = 1 - sW;
    else
        s_hydrocarbon = 1;
    end
    [sG, rs, rv] = calculateHydrocarbonsFromStatusBO(model, st, s_hydrocarbon, x, rs, rv, p);

    
    % We will solve for pressure, water and gas saturation (oil saturation
    % follows via the definition of saturations) and well rates + bhp.
    if model.water
        primaryVars = {'pressure', 'sW', gvar, wellVarNames{:}};
        sat = {sW, 1-sW-sG, sG};
    else
        primaryVars = {'pressure', gvar, wellVarNames{:}};
        sat = {1-sG, sG};
    end
    
    wellSol = DynamicState(wellSol, [wellVarNames, 'wellmap'], [wellVars, wellMap]);
    
    [fp, fpname] = model.FlowPropertyFunctions.getPropertyContainer();
    
    state = DynamicState(state, {'pressure', 's', 'rv', 'rs', 'wellSol', fpname},...
                                {p, sat, rv, rs, wellSol, fp});

%     props = {'PhaseDensities', 'PhaseViscosities', ...
%              'PhaseMobility', 'PhaseComposition', 'PhasePressures'};
    
%     state.evaluatedProperties.Flow = FlowProperties();
end