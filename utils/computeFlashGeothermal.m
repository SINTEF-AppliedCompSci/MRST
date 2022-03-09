function state = computeFlashGeothermal(model, state, varargin)
%Flash geothermal state so that it is in thermodynamic equilibrium

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

    if model.getNumberOfPhases == 2
        state = computePhases(model, state, varargin{:});
    else
        state.flag = ones(model.G.cells.num,1);
    end

    switch model.thermalFormulation
        case 'temperature'
            h = computeEnthalpy(model, state);
            state.enthalpy = h;
        case 'enthalpy'
            T = computeTemperature(model, state);
            state.T = T;
    end
    
end

%-------------------------------------------------------------------------%
function h = computeEnthalpy(model, state)
% Compute enthalpy from a given pressure and temperature
    % Get pressure, internal energy and density
    [p, u, rho] = model.getProps(state, 'pressure'           , ...
                                        'PhaseInternalEnergy', ...
                                        'Density'            );
    % Compute enthalpy
    h = u{1} + p./rho{1};
end

%-------------------------------------------------------------------------%
function T = computeTemperature(model, state)
% Compute temperature for a given pressure and enthalpy
    if isfield(model.fluid, 'T')
        [p, h] = model.getProps(state, 'pressure', 'enthalpy');
        T = model.fluid.T(p, h);
        return
    end
    % Sample variable for Jacobian calculations
    x = model.getProps(state, 'pressure');
    % Set tolerance and maximum iterations
    [tol, itmax] = deal(1e-10, 100);
    % Get current values withou derivatives
    st = value(state);
    T  = model.getProps(st, 'temperature');
    % Solve for pressure
    for i = 1:itmax
        % Set current temperature estimate
        st.T = initVariablesADI(T);
        % Make sure we do not have cahced properties
        st   = model.initStateFunctionContainers(st);
        % Compute residual equation
        res = getTemperatureEquation(model, st);
        if ~isa(res, 'ADI')
            % Temperature can be comuted directly. Ensure that it is
            % consistent, and compute with correct derivatives
            assert(norm(eq, inf) == 0);
            [p, u] = model.getProps(state, 'pressure'           , ...
                                           'PhaseInternalEnergy');
            T = u{1}./model.fluid.CpW(p,T);
            return;
        end
        % Check convergence
        val = norm(value(res), inf)/mean(T);
        if val < tol, break; end
        % Update temperature
        dT = -res.jac{1}\res.val;
        T  = value(T) + dT;
    end
    if ~isa(x, 'ADI'), return; end
    % Compute derivatives using chain rule
    state.T = T;
    % Residual equation with correct Jacobians
    resAD = getTemperatureEquation(model, state);
    % Initialize T as AD
    T      = model.AutoDiffBackend.convertToAD(value(T), x);
    isDiag = isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend');
    d      = full(diag(res.jac{1}));
    for i = 1:numel(resAD.jac)
        if isDiag
            if T.jac{i}.isZero
                T.jac{i} = T.jac{i}.expandZero();
            end
            T.jac{i}.diagonal = -resAD.jac{i}.diagonal./d;
            if all(T.jac{i}.subset == 0)
                T.jac{i}.subset = [];
            end
        else
            T.jac{i} = -res.jac{1}\resAD.jac{i};
        end
    end
end

%-------------------------------------------------------------------------%
function res = getTemperatureEquation(model, state)
    [p, T, u] = model.getProps(state, 'pressure'           , ...
                                      'temperature'        , ...
                                      'PhaseInternalEnergy');
    res = T - u{1}./model.fluid.CpW(p, T);
end

%-------------------------------------------------------------------------%
function state = computePhases(model, state, state0)
    
    [p, h] = model.getProps(state, 'pressure', 'enthalpy');
    
    hL    = model.fluid.hW(p);
    hV    = model.fluid.hG(p);
    pcrit = model.fluid.pcrit;
    hcrit = model.fluid.hcrit;
    
    [isL, isV, isLV] = getPhases(p, h, hL, hV, pcrit, hcrit);
    if nargin > 2 && 0
        [state, h, isL, isV, isLV] = limitState(model, state, state0, h, hL, hV, isL, isV, isLV);
    end
    
    if any(isLV)
        rhoL = model.fluid.rhoW(p, h);
        rhoV = model.fluid.rhoG(p, h);
        sL = rhoV.*(hV - h)./(h.*(rhoL - rhoV) - (hL.*rhoL - hV.*rhoV));
    else
        sL = zeros(numel(value(state.pressure)),1);
    end
    if ~isa(sL, 'ADI') && isa(state.pressure, 'ADI')
        sL = model.AutoDiffBackend.convertToAD(sL, p);
    end

    minSat = 0;
    sL = min(max(sL,minSat),1-minSat);
    sL = sL.*isLV + isL;
    sV = 1 - sL;
    
    flag = 1*isL + 2*isV + 3*isLV;

    state      = model.setProp(state, 'sW', sL);
    state      = model.setProp(state, 'sG', sV);
    state.flag = flag;
end

%-------------------------------------------------------------------------%
function [isL, isV, isLV] = getPhases(p, h, hL, hV, pcrit, hcrit)
    isSC = value(p) > pcrit & value(h) > hcrit;
    isL  = value(h) < hL & ~isSC; % Liquid
    isV  = value(h) > hV & ~isSC; % Vapor
    isV  = isV | isSC;            % Label supercritical as vapor
    isLV = ~isL & ~isV;           % Two-phase
end

%-------------------------------------------------------------------------%
function [state, h, isL, isV, isLV] = limitState(model, state, state0, h, hL, hV, isL, isV, isLV)
    isL0     = state0.flag == 1;
    isV0     = state0.flag == 2;
    isLV0    = state0.flag == 3;
    to2ph    = isLV  & ~isLV0;
    from2ph  = ~isLV & isLV0;
    switched = to2ph | from2ph;
    dh = median(value(h))*0.01;
    if any(switched)
        h0 = h;
        % Case 1: from liquid to two-phase
        l2lv    = value(h) > value(hL) & to2ph;
        h(l2lv) = min(h(l2lv), hL(l2lv) + dh);
        % Case 2: from two-phase to liquid
        lv2l    = value(h) < value(hL) & from2ph;
        h(lv2l) = max(h(lv2l), hL(lv2l) - dh);
        % Case 3: from vapor to two-phase
        v2lv    = value(h) < value(hV) & to2ph & ~l2lv;
        h(v2lv) = max(h(v2lv), hV(v2lv) - dh);
        % Case 4: from two-phase to vapor
        lv2v    = value(h) > value(hV) & from2ph;
        h(lv2v) = min(h(lv2v), hV(lv2v) + dh);
        % Freeze flip-flopping cells
        freeze = state.switchCount > 5;
        isLV(freeze) = isLV0(freeze);
        isL(freeze)  = isL0(freeze);
        isV(freeze)  = isV0(freeze);
        % Set enthalpy
        st = state;
        state = model.setProp(state, 'enthalpy', h);
        state.switchCount = state.switchCount + (switched & ~freeze);
        if 0
            close all
            plot(hL, state.pressure, 'k.');
            hold on
            plot(hV, state.pressure, 'k.');
            plot(state.enthalpy(l2lv), state.pressure(l2lv), '.r')
            plot(state.enthalpy(lv2l), state.pressure(lv2l), '.g')
            plot(state.enthalpy(v2lv), state.pressure(v2lv), '.b')
            plot(state.enthalpy(lv2v), state.pressure(lv2v), '.m')
            plot(st.enthalpy(l2lv), st.pressure(l2lv), 'or')
            plot(st.enthalpy(lv2l), st.pressure(lv2l), 'og')
            plot(st.enthalpy(v2lv), st.pressure(v2lv), 'ob')
            plot(st.enthalpy(lv2v), st.pressure(lv2v), 'om')
            state;
        end
    end
end