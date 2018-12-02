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
        % status 1 oil, no gas  : x = rs, sg = 0    , rv = rvMax
        % status 2 gas, no oil  : x = rv, sg = 1-sw , rs = rsMax
        % status 3 oil and gas  : x = sg, rs = rsMax, rv = rvMax

%     
%     if useAD && false
%         stol = 1e-6;
%         bad_water = double(sW) < stol;% & st{1};
%         if any(bad_water)
% %             sW(bad_water) = stol;
%             sW.val(bad_water) = stol;
%         end
%     end
    
    if isempty(sW)
        sW = 0;
    end
    sG = st{2}.*(1-sW) + st{3}.*x;
    sO = st{1}.*(1-sW) + ~st{1}.*(1 - sW - sG);
%     
%     sO = 1 - sW - sG;

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
    if model.disgas
        rsSat = model.FlowPropertyFunctions.getProperty(model, state, 'RsMax');
        rs = ~st{1}.*rsSat + st{1}.*x;
        rs = rs.*(double(sO) > 0);
        state.rs = rs;
    end
    
    if model.vapoil
        rvSat = model.FlowPropertyFunctions.getProperty(model, state, 'RvMax');
        rv = ~st{2}.*rvSat + st{2}.*x;
        rv = rv.*(double(sG) > 0);
        state.rv = rv;
    end
end