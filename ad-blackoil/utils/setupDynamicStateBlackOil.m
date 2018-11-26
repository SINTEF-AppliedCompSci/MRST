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
    
    if useAD && false
        stol = 1e-6;
        bad_water = double(sW) < stol;% & st{1};
        if any(bad_water)
%             sW(bad_water) = stol;
            sW.val(bad_water) = stol;
        end
    end
    
    if isempty(sW)
        sW = 0;
    end
    
    [sG, rs, rv] = calculateHydrocarbonsFromStatusBO(model, st, 1 - sW, x, rs, rv, p);

    if model.water
        sO = 1 - sW - sG;
    else
        sO = 1 - sG;
    end
    
    if model.vapoil
        % No rv, no so -> zero on diagonal in matrix
        bad_oil = double(sO) == 0 & double(rv) == 0;
        if any(bad_oil)
            sO(bad_oil) = 1 - sW(bad_oil) - double(sG(bad_oil));
        end
    end


    
    % We will solve for pressure, water and gas saturation (oil saturation
    % follows via the definition of saturations) and well rates + bhp.
    if model.water
        primaryVars = {'pressure', 'sW', gvar, wellVarNames{:}};
        sat = {sW, sO, sG};
    else
        primaryVars = {'pressure', gvar, wellVarNames{:}};
        sat = {sO, sG};
    end
    
    wellSol = DynamicState(wellSol, [wellVarNames, 'wellmap'], [wellVars, wellMap]);
    
    [fp, fpname] = model.FlowPropertyFunctions.getPropertyContainer();
    
    state = DynamicState(state, {'pressure', 's', 'rv', 'rs', 'wellSol', fpname},...
                                {p, sat, rv, rs, wellSol, fp});

%     props = {'PhaseDensities', 'PhaseViscosities', ...
%              'PhaseMobility', 'PhaseComposition', 'PhasePressures'};
    
%     state.evaluatedProperties.Flow = FlowProperties();
end