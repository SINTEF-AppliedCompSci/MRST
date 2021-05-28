function props = convertRSRVtoKValue(model, varargin)
% Convert RS/RV b lackoil to K-value. Experimental!
% 
% NOTE:
%   Does not yet support RV.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('np', 100, 'nz', 500, ...
                 'nz_sample', [], 'np_sample', [],...
                 'interpolation', 'linear', ...
                 'minPressure', 300*barsa, ...
                 'maxPressure', 500*barsa);
    opt = merge_options(opt, varargin{:});
    fluid = model.fluid;
    np = 50;
    nz = 50;
    
    np = opt.np;
    nz = opt.nz;
        
    if isempty(opt.np_sample)
        np_s = np;
    else
        np_s = np;
    end
    if isempty(opt.nz_sample)
        nz_s = nz;
    else
        nz_s = nz;
    end
    
    p0 = linspace(opt.minPressure, opt.maxPressure, np)';
    
    if model.disgas
        rs0 = fluid.rsSat(p0);
        rs0 = linspace(1, 2*max(rs0), nz);
    else
        rs0 = linspace(1, 100, nz);
    end


    [P, RS] = meshgrid(p0, rs0);
    SG = 0*RS;
    for pNo = 1:numel(p0)
        if model.disgas
            rssat = fluid.rsSat(P(1, pNo));
        else
            rssat = 0;
        end
        act = RS(:, pNo) > rssat;
        RS(act, pNo) = rssat;
        SG(act, pNo) = linspace(0, 1, nnz(act));
    end
    
    p = P(:);
    sG = SG(:);
    rs = RS(:);

    sO = 1 - sG;
    % Evaluate properties
    [xo, xg, yo, yg, zo, zg, rhoO, rhoG, muO, muG, isSat] = blackOilToMassFraction(model, p, sO, sG, rs, nan);

    
    
    % Equilibrium constants
    Ko = yo./xo;
    Kg = yg./xg;
    
    k_values = struct();
    z0 = linspace(0, 1, nz_s);
    p0 = linspace(p0(1), p0(end), np_s);
    
    [p_s, zg_s] = meshgrid(p0, z0);
    k_values.K_g = interpolateToMesh2D(p, zg, Kg, p_s, zg_s, opt.interpolation);
    k_values.K_o = interpolateToMesh2D(p, zg, Ko, p_s, zg_s, opt.interpolation);
    k_values.freeGas = interpolateToMesh2D(p, zg, double(isSat), p_s, zg_s, 'natural');%, 'nearest');
    
    
    maxK = 1e8;
    k_values.K_g(~isfinite(k_values.K_g)) = maxK;
    k_values.K_o(~isfinite(k_values.K_o)) = maxK;
    
    k_values.z_g = zg_s;
    k_values.z_o = 1 - zg_s;
    k_values.pressure = p_s;
    
    
    % Densities
    [p_s, z_s] = meshgrid(p0, z0);
    properties = struct();
    properties.rhoO = interpolateToMesh2D(p, xo, rhoO, p_s, z_s, opt.interpolation);
    properties.rhoG = interpolateToMesh2D(p, yg, rhoG, p_s, z_s, opt.interpolation);
    
    properties.muO = interpolateToMesh2D(p, xo, muO, p_s, z_s, opt.interpolation);
    properties.muG = interpolateToMesh2D(p, yg, muG, p_s, z_s, opt.interpolation);
    
    properties.pressure = p_s;
    properties.x_g = 1 - z_s;
    properties.x_o = z_s;
    properties.y_g = z_s;
    properties.y_o = 1 - z_s;


    % Viscosities
    
    % Store output
    % K_o(p, z), K_g(p, z)
    % rhoO(p, x), rhoG(p, y)
    % muO(p, x), muG(p, y)

    props = struct();
    props.k_values = k_values;
    props.properties = properties;
    
    
end


function z = interpolateToMesh2D(X, Y, Z, sx, sy, type)
    if nargin < 6
        type = 'linear';
    end
    if isscalar(X)
        [Y, subs] = unique(Y);
    elseif isscalar(Y)
        [X, subs] = unique(X);
    else
        [xy, subs] = unique([X, Y], 'rows');
        X = xy(:, 1);
        Y = xy(:, 2);
    end
    Z = Z(subs);
    
    if max(Y) == min(Y)
        z = interp1(X, Z, sx, type, 'extrap');
    elseif max(X) == min(X)
        z = interp1(Y, Z, sy, type, 'extrap');
    else
        try
            T = scatteredInterpolant(X, Y, Z, type, 'nearest');
        catch
            T = scatteredInterpolant(X, Y, Z, 'linear', 'nearest');
        end
        z = T(sx, sy);
    end
end
















