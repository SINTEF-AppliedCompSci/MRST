function force = getForceProperties(force, model, state)
%Compute properties at the boundary

    is_bc = isfield(force, 'face');
    propsRes         = getPropertiesReservoirSide(model, state, force, is_bc);
    propsForce       = getPropertiesForceSide(model, propsRes, force, is_bc);
    force.propsRes   = propsRes;
    force.propsForce = propsForce;
    force.sat = value(force.propsForce.s);
    % Hack to get correct relperm dimensions when nph = 1
    if is_bc, nc = numel(force.face); else, nc = numel(force.cell); end
    model.G.cells.num = nc;
    [mob, rho] = model.getProps(propsForce, 'Mobility', 'Density');
    force.mob = value(mob);
    force.rho = value(rho);
    
end

%-------------------------------------------------------------------------%
function props  = getPropertiesReservoirSide(model, state, force, is_bc)
    [h, p, T, s, rho, mob] = model.getProps(state, 'enthalpy'                 , ...
                                                   'pressure'                 , ...
                                                   'Temperature'              , ...
                                                   's'                        , ...
                                                   'Density'                  , ...
                                                   'Mobility'                 , ...
                                                   'FluidHeatTransmissibility', ...
                                                   'RockHeatTransmissibility' );
                                               
    if is_bc
        faces = force.face;
        cells = sum(model.G.faces.neighbors(faces,:),2);
    else
        cells = force.cell;
        fpos  = mcolon(model.G.cells.facePos(cells), model.G.cells.facePos(cells+1)-1);
        faces = model.G.cells.faces(fpos);
        nf    = diff(model.G.cells.facePos); nf = nf(cells);
        [ii, jj] = blockDiagIndex(ones(numel(nf),1), nf);
        S        = sparse(ii, jj, 1);
        map      = @(x) (S*x)./nf;
    end
    % Extract BC cell values
    getBCVal = @(v) cellfun(@(v) v(cells), v, 'UniformOutput', false);
    p   = p(cells);
    T   = T(cells);
    h   = h(cells);
    flag = state.flag(cells);
    rho = getBCVal(rho);
    mob = getBCVal(mob);
    % Extract BC face values
    
    if model.dynamicFlowTrans
        force = model.FlowDiscretization.getStateFunction('Transmissibility');
        state = model.FlowDiscretization.evaluateDependencies(model, state, force.dependencies);
        Tf = force.evaluateOnDomain(model, state, true);
        Tf = Tf(faces);
    else
        Tf = model.operators.T_all(faces);
    end
    if model.dynamicHeatTransFluid
        force = model.FlowDiscretization.getStateFunction('FluidHeatTransmissibility');
        state = model.FlowDiscretization.evaluateDependencies(model, state, force.dependencies);
        Thf = force.evaluateOnDomain(model, state, true);
        Thf = Thf(faces);
    else
        Thf = model.operators.Thf_all(faces);
    end
    if model.dynamicHeatTransRock
        force = model.FlowDiscretization.getStateFunction('RockHeatTransmissibility');
        state = model.FlowDiscretization.evaluateDependencies(model, state, force.dependencies);
        Thr = force.evaluateOnDomain(model, state, true);
        Thr = Thr(faces);
    else   
        Thr = model.operators.Thr_all(faces);
    end
    if ~is_bc, [Tf, Thf, Thr] = deal(map(Tf), map(Thf), map(Thr)); end
    if iscell(s)
        s = getBCVal(s);
        s = {s};
    else
        s = s(cells,:);
    end
    props = struct('pressure', p    , ...
                   'enthalpy', h    , ...
                   'T'       , T    , ...
                   's'       , s    , ...
                   'flag'    , flag , ...
                   'rho'     , {rho}, ...
                   'mob'     , {mob}, ...
                   'Tf'      , Tf   , ...
                   'Thr'     , Thr  , ...
                   'Thf'     , Thf  );
end

%-------------------------------------------------------------------------%
function props = getPropertiesForceSide(model, propsRes, force, is_bc) 
    p = getBoundaryPressure(model, propsRes, force, is_bc);
    T = getBoundaryTemperature(model, propsRes, force);
    X = getBoundaryMassFraction(model, force);
    s = force.sat;
    props = struct('pressure', p, ...
                   'X'       , X, ...
                   'T'       , T, ...
                   's'       , s);
    props.components = force.components;
    % Set enthalpy from temperature
    if model.getNumberOfPhases() == 2
        h = model.fluid.h(p,T);
    else
        rho = model.fluid.rhoW(p, T);
        u   = model.fluid.uW(p, T);
        h   = u + p./rho;
    end
    props.enthalpy = h;
    props = computeFlashGeothermal(model, props);
end

%-------------------------------------------------------------------------%
function p = getBoundaryPressure(model, propsRes, force, is_bc)
%     p = model.AutoDiffBackend.convertToAD(zeros(numel(bc.face), 1), propsRes.pressure);
    p = propsRes.pressure;
    if is_bc
        is_pressure = reshape(strcmpi(force.type, 'pressure'), [], 1);
        is_flux     = reshape(strcmpi(force.type, 'flux'), [], 1);
        cells       = sum(model.G.faces.neighbors(force.face,:), 2);
        value       = force.value;
    else
        is_pressure = false;
        is_flux     = true(numel(force.cell), 1);
        cells       = force.cell;
        value       = force.rate;
    end
    % Get pressure at boundaries open to flow
    if any(is_pressure)
        p(is_pressure) = force.value(is_pressure);
    end
    
    if any(is_flux)
        % Flux BCs requires reconstruction. Assuming simple
        G = model.G;
        q = -value(is_flux);
        dzbc = 0;
        if is_bc
            if any(strcmpi(G.type,'topSurfaceGrid'))
                dzbc = model.gravity(3)*(G.cells.z(cells) - G.faces.z(force.face));
            else
                g    = model.getGravityVector();
                dz   = G.faces.centroids(force.face,:) - G.cells.centroids(cells,:);
                dzbc = -dz*g';
            end
            dzbc = dzbc(is_flux);
        end
        [totMob, rhoMob] = deal(0);
        nph = model.getNumberOfPhases();
        for i = 1:nph
            mobi = propsRes.mob{i}(is_flux);
            rhoi = propsRes.rho{i}(is_flux);
            totMob = totMob + mobi;
            rhoMob = rhoMob + mobi.*rhoi;
        end
        T  = propsRes.Tf(is_flux);
        dp = (-q./T + rhoMob.*dzbc)./totMob;
        if ~isa(p, 'ADI') && isa(dp, 'ADI')
            p = model.AutoDiffBackend.convertToAD(p, dp);
        end
        p(is_flux) = p(is_flux) + dp;
    end
end

%-------------------------------------------------------------------------%
function Tbc = getBoundaryTemperature(model, propsRes, force)
    Tbc = model.AutoDiffBackend.convertToAD(force.T, propsRes.pressure);
    is_Hflux = ~isnan(force.Hflux);
    if any(is_Hflux)
        ThR = propsRes.Thr(is_Hflux);
        ThF = propsRes.Thf(is_Hflux);
        Q   = force.Hflux(is_Hflux);
        T   = propsRes.T(is_Hflux);
        dT  = Q./(ThR + ThF);
        Tbc(is_Hflux) = T + dT;
    end
end

%-------------------------------------------------------------------------%
function Xbc = getBoundaryMassFraction(model, bc)
    cnames = model.getComponentNames();
    ix = strcmpi(cnames, 'NaCl');
    Xbc = model.getMassFraction(bc.components);
    if any(ix), Xbc = Xbc(:,ix); else, Xbc = []; end
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}