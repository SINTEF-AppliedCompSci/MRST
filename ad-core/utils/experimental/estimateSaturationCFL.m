function cfl = estimateSaturationCFL(model, state, dt, varargin)
    opt = struct('forces', []);
    opt = merge_options(opt, varargin{:});
    
    pv = model.operators.pv;
    
    [F, Q] = getFlowValue(model, state);
    
    v = state.flux(model.operators.internalConn, :);
    vT = sum(v, 2);
    
    xflow = ~(all(v >= 0, 2) | all(v <= 0, 2));
    
    l = model.operators.N(:, 1);
    r = model.operators.N(:, 2);
    
    flag = vT > 0;
    df_face = upstream(model, F, flag, xflow, l, r);
    
    rate_face = abs(vT).*abs(df_face);
    if ~isempty(Q)
        T = model.operators.T;
        cap_face = upstream(model, Q, flag, xflow, l, r);
        rate_face = rate_face + 2.*T.*cap_face;
    end
    
    nc = model.G.cells.num;
    % Accumulate into cell if flow is outgoing, or we have any kind of
    % cross-flow.
    rate_cell = accumarray(l, rate_face.*( flag | xflow), [nc, 1]) +...
                accumarray(r, rate_face.*(~flag | xflow), [nc, 1]);
    
    cfl = (dt./pv).*rate_cell;
end

function df_face = upstream(model, F, flag, xflow, l, r)
    df_face = model.operators.faceUpstr(flag, F);
    df_face(xflow) = max(abs(F(l(xflow))), abs(F(r(xflow))));
end

function [F, Q] = getFlowValue(model, state)
    nph = model.water + model.oil + model.gas;
    nc = model.G.cells.num;
    
    if nph == 1
        F = ones(nc, 1);
    elseif nph == 2
        if model.oil
            if model.gas
                [F, Q] = getFlowMagnitudeOG(model, state);
            elseif model.water
                [F, Q] = getFlowMagnitudeOW(model, state);
            end
        else
            error('Not supported');
        end
    else
        [F, Q] = getFlowMagnitudeWOG(model, state);
    end
    
    F = findLargestEigenvalue(F);
    Q = findLargestEigenvalue(Q);
end

function [F, Q] = getFlowMagnitudeWOG(model, state)
    assert(model.water && model.oil && model.gas)
    
    [sW, sG] = model.getProps(state, 'sW', 'sG');
    [sW, sG] = initVariablesADI(sW, sG);
    sO = 1 - sW - sG;
    
    [krW, krO, krG] = model.evaluateRelPerm({sW, sO, sG});
    muW = getViscosity(model, state, 'w');
    muO = getViscosity(model, state, 'o');
    muG = getViscosity(model, state, 'g');
    
    mobW = krW./muW;
    mobO = krO./muO;
    mobG = krG./muG;
    
    mobT = mobO + mobG + mobW;
    
    f_g = mobG./mobT;
    f_w = mobW./mobT;
    
    J = @(x, i) full(diag(x.jac{i}));
    F = [J(f_w, 1), J(f_g, 1), J(f_w, 2), J(f_g, 2)];
    
    f = model.fluid;
    hasOW = isfield(f, 'pcOW');
    hasOG = isfield(f, 'pcOG');
    
    if hasOW || hasOG
    	Q = zeros(model.G.cells.num, 4);
        mobOW = double(mobO.*mobW)./double(mobT);
        mobOG = double(mobO.*mobG)./double(mobT);
        mobWG = double(mobW.*mobG)./double(mobT);
        if hasOW
            pcow = f.pcOW(sW);
            dow = J(pcow, 1);
            Q(:, 1) =  dow.*mobWG;
            Q(:, 2) = -dow.*(mobWG + mobOW);
        end
        if hasOG
            pcog = f.pcOG(sG);
            dog = J(pcog, 2);
            Q(:, 3) =  dog.*(mobOG + mobWG);
            Q(:, 4) = -dog.*mobWG;
        end
    else
        Q = [];
    end
end



function [F, Q] = getFlowMagnitudeOG(model, state)
    assert(~model.water && model.oil && model.gas)
    
    sG = model.getProp(state, 'sG');
    sG = initVariablesADI(sG);
    sO = 1 - sG;
    
    [krO, krG] = model.evaluateRelPerm({sO, sG});
    muG = getViscosity(model, state, 'g');
    muO = getViscosity(model, state, 'o');

    mobG = krG./muG;
    mobO = krO./muO;
    
    mobT = mobO + mobG;
    
    f_g = mobG./mobT;
    J = @(x) full(diag(x.jac{1}));

    F = J(f_g);
    
    
    hasOG = isfield(f, 'pcOG');
    if hasOG
    	Q = zeros(model.G.cells.num, 1);
        mobOG = double(mobO.*mobG)./double(mobT);
        if hasOG
            pcog = f.pcOG(sG);
            dog = J(pcog);
            Q(:, 1) = dog.*mobOG;
        end
    else
        Q = [];
    end

end

function [F, Q] = getFlowMagnitudeOW(model, state)
    assert(model.water && model.oil && ~model.gas)
    
    sW = model.getProp(state, 'sW');
    sW = initVariablesADI(sW);
    sO = 1 - sW;
    
    [krW, krO] = model.evaluateRelPerm({sW, sO});
    muW = getViscosity(model, state, 'w');
    muO = getViscosity(model, state, 'o');

    mobW = krW./muW;
    mobO = krO./muO;
    
    mobT = mobO + mobW;
    
    f_w = mobW./mobT;
    
    J = @(x) full(diag(x.jac{1}));
    F = J(f_w);
    
    hasOW = isfield(f, 'pcOW');
    
    if hasOW 
    	Q = zeros(model.G.cells.num, 1);
        mobOW = double(mobO.*mobW)./double(mobT);
        pcow = f.pcOW(sW);
        dow = J(pcow);
        Q(:, 1) =  dow.*mobOW;
    else
        Q = [];
    end

end


function mu = getViscosity(model, state, name)
    isOil = strcmpi(name, 'o');
    isGas = strcmpi(name, 'g');
    
    p = model.getProps(state, 'pressure');
    if isa(model, 'ThreePhaseCompositionalModel') && (isOil || isGas)
        pm = model.EOSModel.PropertyModel;
        isLiquid = strcmpi(name, 'o');
        if isLiquid
            [xy, T, Z] = model.getProps(state, 'x', 'T', 'Z_L');
        else
            [xy, T, Z] = model.getProps(state, 'y', 'T', 'Z_V');
        end
        mu = computeViscosity(pm, p, xy, Z, T, isLiquid);
    else
        fn_name = ['mu', upper(name)];
        mu = model.fluid.(fn_name)(p);
    end

end


function F = findLargestEigenvalue(F)
    if size(F, 2) > 1
        nc = size(F, 1);
        nph = sqrt(size(F, 2)) + 1;
        tmp = zeros(nc, 1);
        for i = 1:nc
            M = reshape(F(i, :), nph-1, nph-1);
            tmp(i) = max(abs(eig(M)));
        end
        F = tmp;
    else
        F = abs(F);
    end
end