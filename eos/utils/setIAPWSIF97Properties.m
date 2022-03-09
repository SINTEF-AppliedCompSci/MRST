function fluid = setIAPWSIF97Properties(fluid, varargin)
%Set IAPWSIF97 water properties to fluid

    opt = struct('prange', [0.1, 100]*mega*Pascal  , ...
                 'Trange', [275.16, 1073.15]*Kelvin);
    [opt, extra] = merge_options(opt, varargin{:});
    % Set pressure, temperature and enthalpy range
    [prange, Trange, hrange] = getVariableRange(opt);
    % Set critical pressure and temperature
    fluid = setCriticalPoint(fluid);
    % Get input data to generate Coolprop table
    in = getInputData(fluid, prange, Trange, hrange);
    % Get tables and postprocess
    tab = struct();
    n   = size(in,1);
    for i = 1:n
        name       = in{i,1};
        tab.(name) = generateCoolPropTable(in{i,2:end}, extra{:}, 'save', true);
        tab.(name) = postprocessPVTTable(tab.(name), name);
    end
    % Set function handles to fluid from tables
    for i = 1:n
        name         = in{i,1};
        fluid.(name) = getPVTFunction(tab.(name), name);
    end
    % Set surface densities (may be used some places in the code)
    fluid = setSurfaceDensities(fluid);
    % Overwrite thermal conductivity. TODO: Fix this table
    Watt = joule/second;
    fluid.lambdaF = @(varargin) varargin{1}*0 + 0.6*Watt/(meter*Kelvin);
    
end

%-------------------------------------------------------------------------%
function phases = getPhases(fluid)
    names  = fieldnames(fluid)';
    assert(~any(strcmpi(names, 'krO')), 'Fluid cannot have an oil phase');
    water  = any(strcmpi(names, 'krW'));
    gas    = any(strcmpi(names, 'krG'));
    phases = 'WG';
    phases = phases([water, gas]);
end

%-------------------------------------------------------------------------%
function [p, T, h] = getVariableRange(opt)
    % Maximum valid range
    p = [  0.1 , 100    ]*mega*Pascal;
    T = [275.16, 1073.15];
    % Check pressure range
    if isempty(opt.prange), opt.prange = p; end
    if opt.prange(1) < min(p(1)) || opt.prange(2) > min(p(2))
        warning('Pressure outside valid range. Restricting to valid range.')
        p = [max(p(1), opt.prange(1)), min(p(2), opt.prange(2))];
    end
    % Check temperature range
    if isempty(opt.Trange), opt.Trange = T; end
    if opt.prange(1) < min(p(1)) || opt.prange(2) > min(p(2))
        warning('Pressure outside valid range. Restricting to valid range.')
        T = [max(T(1), opt.Trange(1)), min(T(2), opt.Trange(2))];
    end
    % Set enthalpy range
    h = [py.CoolProp.CoolProp.PropsSI('H', 'P', p(1), 'T', T(1), 'IF97::Water'), ...
         py.CoolProp.CoolProp.PropsSI('H', 'P', p(2), 'T', T(2), 'IF97::Water')];
end

%-------------------------------------------------------------------------%
function fluid = setCriticalPoint(fluid)
    fluid.pcrit = py.CoolProp.CoolProp.PropsSI('Pcrit', 'IF97::Water');    
    fluid.Tcrit = py.CoolProp.CoolProp.PropsSI('Tcrit', 'IF97::Water');  
    fluid.hcrit = py.CoolProp.CoolProp.PropsSI('H', 'P', fluid.pcrit, 'T', fluid.Tcrit, 'IF97::Water');
end

%-------------------------------------------------------------------------%
function in = getInputData(fluid, prange, Trange, hrange)
    % Get active phases
    phases = getPhases(fluid);
    nph    = numel(phases);
    switch nph
        case 1
            name = @(prefix) [prefix, phases];
            in = { ...
                % fluid.n    % Coolprop name   % Input 1   % Input 2      % Range
                'T'        , 'temperature'   , 'pressure', 'enthalpy'   , prange, hrange;
                name('rho'), 'density'       , 'pressure', 'enthalpy'   , prange, hrange;
                name('mu') , 'viscosity'     , 'pressure', 'temperature', prange, Trange;
                name('h')  , 'enthalpy'      , 'pressure', 'temperature', prange, Trange;
                name('u')  , 'internalenergy', 'pressure', 'enthalpy'   , prange, hrange;
                name('Cp') , 'heatcapacity'  , 'pressure', 'temperature', prange, [0,0] ;
                'lambdaF'  , 'conductivity'  , 'pressure', 'temperature', prange, Trange;
            };
        case 2
            % Input data
            in = { ...
                % fluid.n   % Coolprop name   % Input 1   % Input 2      % Range
                'T'       , 'temperature'   , 'pressure', 'enthalpy'   , prange, hrange;
                'Tboiling', 'temperature'   , 'pressure', 'quality'    , prange, [0,0] ;
                'rho'     , 'density'       , 'pressure', 'enthalpy'   , prange, hrange;
                'rhoW'    , 'density'       , 'pressure', 'quality'    , prange, [0,0] ;
                'rhoG'    , 'density'       , 'pressure', 'quality'    , prange, [1,1] ;
                'mu'      , 'viscosity'     , 'pressure', 'enthalpy'   , prange, hrange;
                'muW'     , 'viscosity'     , 'pressure', 'quality'    , prange, [0,0] ;
                'muG'     , 'viscosity'     , 'pressure', 'quality'    , prange, [1,1] ;
                'h'       , 'enthalpy'      , 'pressure', 'temperature', prange, Trange;
                'hW'      , 'enthalpy'      , 'pressure', 'quality'    , prange, [0,0] ;
                'hG'      , 'enthalpy'      , 'pressure', 'quality'    , prange, [1,1] ;
                'u'       , 'internalenergy', 'pressure', 'enthalpy'   , prange, hrange;
                'uW'      , 'internalenergy', 'pressure', 'quality'    , prange, [0,0] ;
                'uG'      , 'internalenergy', 'pressure', 'quality'    , prange, [1,1] ;
                'Cp'      , 'heatcapacity'  , 'pressure', 'temperature', prange, Trange;
                'CpW'     , 'heatcapacity'  , 'pressure', 'quality'    , prange, [0,0] ;
                'CpG'     , 'heatcapacity'  , 'pressure', 'quality'    , prange, [1,1] ;
                'lambdaF' , 'conductivity'  , 'pressure', 'temperature', prange, Trange;
                'lambdaFW', 'conductivity'  , 'pressure', 'quality'    , prange, [0,0] ;
                'lambdaFG', 'conductivity'  , 'pressure', 'quality'    , prange, [1,1] ;
            };
    end

end

%-------------------------------------------------------------------------%
function tab = postprocessPVTTable(tab, name)
    args = {};
    if numel(tab.y) > 1
        switch name
            case 'T'
                args = {'fixPatch' , true                     , ...
                        'rangeX'   , [3.6    , 4.8       ]*1e6, ...
                        'rangeY'   , [2.8*1e6, max(tab.y)]    , ...
                        'primaryAx', 'x'                      };
            case 'mu'
                args = {'fixZeros', true};
        end
        tab = postprocessTable(tab, args{:});
    end
end

%-------------------------------------------------------------------------%
function fun = getPVTFunction(tab, name)
    if numel(tab.y) > 1
        fun = getMultiDimInterpolator({tab.x, tab.y}, tab.data, 'nearest');
        fun = @(varargin) fun(varargin{1:2});
    else
        fun = @(varargin) interpTableMEX(tab.x, tab.data, varargin{1});
    end
end

%-------------------------------------------------------------------------%
function fluid = setSurfaceDensities(fluid)
    phases = getPhases(fluid);
    nph    = numel(phases);
    K0 = 273.15 + 20; p0 = 1*atm;
    if nph == 1
        name = @(prefix) [prefix, phases];
    else
        name = @(prefix) prefix;
    end
    h0 = fluid.(name('h'))(p0, K0);
    for i = 1:nph
        fluid.(['rho', phases(i), 'S']) = fluid.(name('rho'))(p0, h0);
    end
end

%-------------------------------------------------------------------------%
function flag = getPhaseFlag(tab, hL, hV, pcrit)
    n = numel(tab.x);
    p = tab.x;
    hL = getPVTFunction(hL); hL = hL(p);
    hV = getPVTFunction(hV); hV = hV(p);
    p = repmat(p, 1, n);
    hV = repmat(hV, 1, n);
    hL = repmat(hL, 1, n);
    h  = repmat(tab.y', n, 1);
    [isL, isV, isLV] = getPhaseFlagGeothermal(p, h, hL, hV, pcrit);
    flag = 1*isL + 2*isV + 3*isLV;
end

%-------------------------------------------------------------------------%
function tab = splitTable(tab, flag)

    close all
    if numel(tab.y) == 1, return; end
    [x, y] = ndgrid(tab.x, tab.y);
    
    
    isL = flag == 1;
    xL = x(isL); yL = y(isL); zL = tab.data(isL);
    F = scatteredInterpolant(xL, yL, zL, 'linear', 'linear');
    zL = F(x,y);
    
    zmin = min(tab.data(:));
    zmax = max(tab.data(:));
    cl = linspace(zmin, zmax, 10);
    figure(), contourf(x, y, tab.data, cl); caxis([zmin, zmax])
    figure(), contourf(x, y, zL, cl); caxis([zmin, zmax])
    figure(), contourf(x, y, isL*1, 1)

    isV = flag == 2;
    xV = x(isV); yV = y(isV); zV = tab.data(isV);
    F = scatteredInterpolant(xV, yV, zV, 'linear', 'linear');
    zV = F(x,y);
    figure(), contourf(x, y, zV, cl); caxis([zmin, zmax])
    figure(), contourf(x, y, isV*1, 1)

    isLV = flag == 3;
    xLV = x(isLV); yLV = y(isLV); zLV = tab.data(isLV);
    F = scatteredInterpolant(xLV, yLV, zLV, 'linear', 'linear');
    zLV = F(x,y);
    figure(), contourf(x, y, zLV, cl); caxis([zmin, zmax])
    figure(), contourf(x, y, isLV*1, 1)

    tab.L  = zL;
    tab.V  = zV;
    tab.LV = zLV;

end

%-------------------------------------------------------------------------%
function tab = smoothBoilingCurve(fluid, tab, in, prange, hrange)
    return
    n = size(in,1);
    np = numel(tab.T.x);
    p = linspace(prange(1), prange(2), np);
    h = linspace(hrange(1), hrange(2), np);
    [p, h] = ndgrid(p, h);
    hL = reshape(fluid.hW(p(:)), np, np);
    hV = reshape(fluid.hG(p(:)), np, np);
    dh = mean(mean(h))*0.05;
    isSC = p > fluid.pcrit;
    reg = (abs(h - hL) < dh | abs(h - hV) < dh) & ~isSC;
    v = ones(np^2, 1);
    A = spdiags([v,v,v,v,v], [-np, 1, 0, 1, np], np^2, np^2);
    reg = reg(:);
    for i = 1:5
        reg = reg | A*reg > 0;
    end
    reg = reshape(reg, np, np);
    
    rhoL = reshape(fluid.rhoW(p), np, np);
    rhoV = reshape(fluid.rhoG(p), np, np);
    hL   = reshape(fluid.hW(p), np, np);
    hV   = reshape(fluid.hG(p), np, np);
    sL   = rhoV.*(hV - h)./(h.*(rhoL - rhoV) - (hL.*rhoL - hV.*rhoV));
    reg = reg & (sL >= 0 & sL <= 1);
    sV   = 1 - sL;
    
    muL = reshape(fluid.muW(p), np, np);
    muV = reshape(fluid.muG(p), np, np);
    
    mu = (sL.*rhoL.*muL + sV.*rhoV.*muV)./(sL.*rhoL + sV.*rhoV);
    tab.mu.data(reg) = muL(reg);
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