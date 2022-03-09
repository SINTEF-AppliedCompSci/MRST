function bc = getBCProperties(bc, model, state)
%Compute properties at the boundary

    propsRes = getReservoirProperties(model, state, bc);
    propsBC  = getBoundaryProperties(model, propsRes, bc);
    bc.propsRes = propsRes;
    bc.propsBC  = propsBC;
    bc.sat = value(bc.propsBC.s);
    % Hack to get correct relperm dimensions when nph = 1
    model.G.cells.num = numel(bc.face);
    [mob, rho] = model.getProps(propsBC, 'Mobility', 'Density');
    bc.mob = value(mob);
    bc.rho = value(rho);
    
end

%-------------------------------------------------------------------------%
function props  = getReservoirProperties(model, state, bc)
    [h, p, T, s, rho, mob] = model.getProps(state, 'enthalpy'                 , ...
                                                   'pressure'                 , ...
                                                   'Temperature'              , ...
                                                   's'                        , ...
                                                   'Density'                  , ...
                                                   'Mobility'                 , ...
                                                   'FluidHeatTransmissibility', ...
                                                   'RockHeatTransmissibility' );
    faces = bc.face;
    cells = sum(model.G.faces.neighbors(faces,:),2);
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
        prop = model.FlowDiscretization.getStateFunction('Transmissibility');
        Tf = prop.evaluateOnDomain(model, state, true);
        Tf = Tf(faces);
    else
        Tf = model.operators.T_all(faces);
    end
    if model.dynamicHeatTransFluid
        prop = model.FlowDiscretization.getStateFunction('FluidHeatTransmissibility');
        Thf = prop.evaluateOnDomain(model, state, true);
        Thf = Thf(faces);
    else
        Thf = model.operators.Thf_all(faces);
    end
    if model.dynamicHeatTransRock
        prop = model.FlowDiscretization.getStateFunction('RockHeatTransmissibility');
        Thr = prop.evaluateOnDomain(model, state, true);
        Thr = Thr(faces);
    else   
        Thr = model.operators.Thr_all(faces);
    end
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
function props = getBoundaryProperties(model, propsRes, bc) 
    p = getBoundaryPressure(model, propsRes, bc);
    T = getBoundaryTemperature(model, propsRes, bc);
    X = getBoundaryMassFraction(model, bc);
    s = zeros(numel(bc.face), model.getNumberOfPhases());
    props = struct('pressure', p, ...
                   'X'       , X, ...
                   'T'       , T, ...
                   's'       , s);
    props.components = bc.components;
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
function p = getBoundaryPressure(model, propsRes, bc)
%     p = model.AutoDiffBackend.convertToAD(zeros(numel(bc.face), 1), propsRes.pressure);
    p = propsRes.pressure;
    is_pressure = reshape(strcmpi(bc.type, 'pressure'), [], 1);
    is_flux     = reshape(strcmpi(bc.type, 'flux'), [], 1);
    % Get pressure at boundaries open to flow
    if any(is_pressure)
        p(is_pressure) = bc.value(is_pressure);
    end
    cells = sum(model.G.faces.neighbors(bc.face,:), 2);
    if any(is_flux)
        % Flux BCs requires reconstruction. Assuming simple
        G = model.G;
        q = -bc.value(is_flux);
        if any(strcmpi(G.type,'topSurfaceGrid'))
            dzbc = model.gravity(3)*(G.cells.z(cells) - G.faces.z(bc.face));
        else
            g    = model.getGravityVector();
            dz   = G.faces.centroids(bc.face,:) - G.cells.centroids(cells,:);
            dzbc = -dz*g';
        end
        dzbc = dzbc(is_flux);
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
        p(is_flux) = p(is_flux) + dp;
    end
end

%-------------------------------------------------------------------------%
function Tbc = getBoundaryTemperature(model, propsRes, bc)
    Tbc = model.AutoDiffBackend.convertToAD(bc.T, propsRes.pressure);
    is_Hflux = ~isnan(bc.Hflux);
    if any(is_Hflux)
        ThR = propsRes.Thr;
        ThF = propsRes.Thf;
        T = propsRes.T;
        dT            = bc.Hflux(is_Hflux)./(ThR(is_Hflux) + ThF(is_Hflux));
        Tbc(is_Hflux) = T(is_Hflux) + dT;
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
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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